function varargout=q4_outputs(action,varargin)
action=string(action);
assert(isscalar(action),'CUMCM:q4:OutputAction','Invalid output action.');
switch action
    case "write_all"
        [varargout{1:nargout}]=write_all(varargin{:});
    case "write_tables"
        varargout{1}=write_q4_outputs(varargin{:});
    case "write_result4"
        varargout{1}=write_q4_result4(varargin{:});
    case "plot_results"
        varargout{1}=plot_q4_results(varargin{:});
    case "export_boundary"
        varargout{1}=export_q4_boundary_deliverables(varargin{:});
    otherwise
        error('CUMCM:q4:OutputAction','Unknown output action: %s',action);
end
end

function paths=write_all(cfg,candidate,summary,in)
paths=write_q4_outputs(cfg,candidate,summary,in);
paths=[paths,plot_q4_results(cfg,candidate,summary,in)];
end

function paths=write_q4_outputs(cfg,candidate,summary,in)
f=candidate.fields;p=q4_compute("parameters",cfg,candidate.parameters);
extra=q4_extra_diagnostics(candidate,in,p);
folder=fullfile(cfg.resultsDir,'q4');if ~isfolder(folder),mkdir(folder);end
paths=strings(1,0);
names=["时间_s",compose("半径_%.1fcm",f.distance_cm),"药材表面"];
put('q4_moisture_unrounded.csv',[f.time_s,f.C],names);
put('q4_temperature_unrounded.csv',[f.time_s,f.temperature_C],names);
put('q4_radius_and_mask.csv',[f.time_s,f.radius_cm,double(f.mask),f.time_s>14400,f.time_s>259200], ...
    ["时间_s","表面半径_cm",compose("域内_%.1fcm",f.distance_cm),"环境外推","半径外推"]);
selected=find(mod(round(f.time_s*1000),21600000)==0);selected=unique([selected;numel(f.time_s)]);
paperDistance=(0:floor(max(f.radius_cm(selected))/.5))*.5;
[tf,ix]=ismember(round(paperDistance*10),round(f.distance_cm*10));assert(all(tf));
put('q4_table6_moisture.csv',[f.time_s(selected)/3600,f.C(selected,[ix,size(f.C,2)])], ...
    ["时间_h",compose("半径_%.1fcm",paperDistance),"药材表面"]);
put('q4_drying_progress.csv',extra.progress,["时间_s","平均含水率_kgkg","达标体积分数","最大含水率_kgkg","环境外推","半径外推"]);
put('q4_mechanism_decomposition.csv',extra.mechanism, ...
    ["时间_s","材料坐标","表面半径_cm","含水率_kgkg","温度_K","参考扩散系数_s_inv", ...
    "总对数变化","含水率贡献","温度贡献","收缩贡献","恒等式闭合残差"]);
event=rmfield(candidate.event,{'stop_state','root_state'});
putjson('q4_event_summary.json',event);
me=candidate.mean_event;
if ~isempty(me),me=rmfield(me,{'stop_state','root_state'});end
putjson('q4_mean_event_summary.json',me);
putjson('q4_validation_summary.json',summary);
save(fullfile(folder,'q4_solution.mat'),'candidate','summary','-v7.3');
paths(end+1)="results/q4/q4_solution.mat";
paths(end+1)=write_q4_result4(cfg,f);
paths=[paths,write_q4_environment_evidence(cfg,in,"results/q4")];
    function put(name,values,headers)
        file=fullfile(folder,name);
        fid=fopen(file,'w','n','UTF-8');assert(fid>0);cleanup=onCleanup(@()fclose(fid));
        fprintf(fid,'%s\n',strjoin(headers,','));
        for row=1:size(values,1)
            cellstrs=strings(1,size(values,2));
            for col=1:size(values,2)
                if isfinite(values(row,col))
                    if string(name)=="q4_table6_moisture.csv"
                        cellstrs(col)=sprintf('%.4f',values(row,col));
                    else
                        cellstrs(col)=sprintf('%.17g',values(row,col));
                    end
                end
            end
            fprintf(fid,'%s\n',strjoin(cellstrs,','));
        end
        paths(end+1)="results/q4/"+name;
    end
    function putjson(name,value)
        fid=fopen(fullfile(folder,name),'w','n','UTF-8');assert(fid>0);
        cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
        paths(end+1)="results/q4/"+name;
    end
end


function path=write_q4_result4(cfg,fields)
source=fullfile(cfg.rawDataDir,'附件3','result4.xlsx');original=hashFile(source);
folder=fullfile(cfg.processedDataDir,'q4');if ~isfolder(folder),mkdir(folder);end
work=fullfile(cfg.projectRoot,'artifacts','q4','output_work');if ~isfolder(work),mkdir(work);end
parts=string(tempname(work));mkdir(parts);copied=string(tempname(folder))+".xlsx";
copyfile(source,copied);unzip(copied,parts);
count=numel(fields.time_s);ncol=size(fields.C,2)+1;
assert(count+1<=1048576&&ncol<=26);
stylesPath=fullfile(parts,'xl','styles.xml');styles=fileread(stylesPath);
token=regexp(styles,'<cellXfs count="(\d+)">(.*?)</cellXfs>','tokens','once');assert(~isempty(token));
xfs=regexp(token{2},'<xf\s[^>]*>.*?</xf>','match');assert(numel(xfs)==str2double(token{1}));
numericStyle=regexprep(xfs{2},'numFmtId="\d+"','numFmtId="4"');
ids=regexp(styles,'numFmtId="(\d+)"','tokens');ids=cellfun(@(x)str2double(x{1}),ids);newId=max([163,ids])+1;
numericStyle=regexprep(numericStyle,'numFmtId="\d+"',sprintf('numFmtId="%d"',newId));
numericStyle=regexprep(numericStyle,' applyNumberFormat="[^"]*"','');
numericStyle=replace(numericStyle,'<xf ','<xf applyNumberFormat="1" ');
fmt=sprintf('<numFmt numFmtId="%d" formatCode="0.0000"/>',newId);
fmts=regexp(styles,'<numFmts count="(\d+)">(.*?)</numFmts>','tokens','once');
if isempty(fmts)
    styles=regexprep(styles,'(<styleSheet[^>]*>)',['$1<numFmts count="1">',fmt,'</numFmts>'],'once');
else
    styles=regexprep(styles,'<numFmts count="\d+">.*?</numFmts>', ...
        sprintf('<numFmts count="%d">%s%s</numFmts>',str2double(fmts{1})+1,fmts{2},fmt),'once');
end
numericIndex=numel(xfs);
styles=regexprep(styles,'<cellXfs count="\d+">.*?</cellXfs>', ...
    sprintf('<cellXfs count="%d">%s%s</cellXfs>',numericIndex+1,token{2},numericStyle),'once');
writeText(stylesPath,styles);
xmlPath=fullfile(parts,'xl','worksheets','sheet1.xml');xml=fileread(xmlPath);
assert(~contains(xml,'<mergeCells')&&~contains(xml,'<f>'));
rows=cell(count+1,1);header=cell(1,ncol);
header{1}='<c r="A1" s="1" t="inlineStr"><is><t>时间&#92;到药材中心的距离</t></is></c>';
for col=2:ncol-1
    header{col}=sprintf('<c r="%s1" s="1"><v>%.1f</v></c>',char(64+col),fields.distance_cm(col-1));
end
header{ncol}=sprintf('<c r="%s1" s="1" t="inlineStr"><is><t>药材表面</t></is></c>',char(64+ncol));
rows{1}=['<row r="1" ht="14.1" customHeight="1">',strjoin(header,''),'</row>'];
values=round(fields.C,4);
for k=1:count
    cells=cell(1,ncol);row=k+1;
    cells{1}=sprintf('<c r="A%d" s="1"><v>%.17g</v></c>',row,fields.time_s(k));
    for col=2:ncol
        if isfinite(values(k,col-1))
            cells{col}=sprintf('<c r="%s%d" s="%d"><v>%.4f</v></c>',char(64+col),row,numericIndex,values(k,col-1));
        else
            cells{col}=sprintf('<c r="%s%d" s="%d"/>',char(64+col),row,numericIndex);
        end
    end
    rows{row}=sprintf('<row r="%d" ht="14.1" customHeight="1">%s</row>',row,strjoin(cells,''));
end
xml=regexprep(xml,'<sheetData>.*?</sheetData>',['<sheetData>',strjoin(rows,''),'</sheetData>']);
xml=regexprep(xml,'<dimension ref="[^"]*"/>',sprintf('<dimension ref="A1:%s%d"/>',char(64+ncol),count+1));
xml=replace(xml,'min="2" max="6"',sprintf('min="2" max="%d"',ncol));
xml=replace(xml,'min="7" max="16384"',sprintf('min="%d" max="16384"',ncol+1));writeText(xmlPath,xml);
candidate=string(tempname(work))+".xlsx";q4_pack_ooxml(parts,candidate);
assert(isequal(sheetnames(candidate),"Sheet1"));
head=readcell(candidate,'Sheet','Sheet1','Range',sprintf('A1:%s1',char(64+ncol)),'UseExcel',false);
assert(string(head{1})=="时间\到药材中心的距离"&&string(head{end})=="药材表面");
matrix=readmatrix(candidate,'Sheet','Sheet1','Range',sprintf('A2:%s%d',char(64+ncol),count+1),'UseExcel',false);
assert(isequal(matrix(:,1),fields.time_s)&&isequaln(matrix(:,2:end),values));
assert(strcmp(hashFile(source),original),'CUMCM:q4:TemplateMutation','Raw template changed.');
destination=fullfile(folder,'result4.xlsx');
if isfile(destination),copyfile(destination,string(tempname(work))+"_previous.xlsx");end
movefile(candidate,destination,'f');path="data/processed/q4/result4.xlsx";
end
function writeText(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',value);
end
function value=hashFile(path)
bytes=java.nio.file.Files.readAllBytes(java.io.File(char(path)).toPath());
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end


function paths=write_q4_environment_evidence(cfg,in,relativeFolder)
folder=fullfile(cfg.projectRoot,relativeFolder);if ~isfolder(folder),mkdir(folder);end
paths=[string(relativeFolder)+"/q4_environment_windows.csv", ...
    string(relativeFolder)+"/q4_environment_selection.json"];
writetable(in.environment_selection.windows,fullfile(cfg.projectRoot,paths(1)),'Encoding','UTF-8');
fid=fopen(fullfile(cfg.projectRoot,paths(2)),'w','n','UTF-8');assert(fid>0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(in.environment_selection.selection,'PrettyPrint',true));
end


function extra=q4_extra_diagnostics(candidate,in,p)
mesh=candidate.mesh;
times=[0;candidate.time_s;candidate.event.root_s;candidate.event.stop_time_s];
states=[repmat([p.T0_K;p.C0],mesh.N+1,1),candidate.states, ...
    candidate.event.root_state,candidate.event.stop_state];
if ~isempty(candidate.mean_event)
    times=[times;candidate.mean_event.root_s;candidate.mean_event.stop_time_s];
    states=[states,candidate.mean_event.root_state,candidate.mean_event.stop_state];
end
[times,ix]=unique(times);states=states(:,ix);progress=zeros(numel(times),6);
mechanism=zeros(3*numel(times),11);
for k=1:numel(times)
    geo=q4_compute("interpolate_radius",times(k),in,p.geometry);bc=q4_compute("interpolate_ambient",times(k),in,p.environment);
    if times(k)==0,bc.moisture_kgkg=p.C0;bc.temperature_K=p.T0_K;end
    d=q4_compute("progress_metrics",states(:,k),mesh,geo,bc,p,true);
    m=q4_compute("spatial_maximum",states(:,k),mesh,geo,bc,p);
    progress(k,:)=[times(k),d.mean_C,d.dry_fraction,m.maximum_C,times(k)>14400,times(k)>259200];
    rows=(k-1)*3+(1:3);
    mechanism(rows,:)=[repmat(times(k),3,1),[0;.5;1],repmat(100*geo.R_m,3,1), ...
        d.C_at_material_points,d.T_at_material_points,d.Lambda,d.total_log, ...
        d.moisture_log,d.temperature_log,d.geometry_log,d.closure_residual];
end
extra=struct('progress',progress,'mechanism',mechanism);
end


function paths=plot_q4_results(cfg,candidate,summary,in)
f=candidate.fields;p=q4_compute("parameters",cfg,candidate.parameters);extra=q4_extra_diagnostics(candidate,in,p);
stems=["raw_q4_radius","result_q4_moisture_field","result_q4_drying_event","result_q4_ablation"];
types=["raw","result","result","result"];paths=strings(1,0);
blue=[0,.447,.698];orange=[.835,.369,0];green=[0,.62,.45];
figureIds=1:4;
if isfield(summary,'primary_only')&&summary.primary_only,figureIds=[2,3];end
for k=figureIds
    fig=figure('Visible','off','Color','w');cleanup=onCleanup(@()close(fig));
    apply_publication_style(fig,'en','report');
    assert(any(string(listfonts)=="Microsoft YaHei UI"),'CUMCM:q4:FontMissing','Chinese figure font unavailable.');
    set(fig,'DefaultAxesFontName','Microsoft YaHei UI','DefaultTextFontName','Microsoft YaHei UI');
    set(fig,'DefaultTextInterpreter','none','DefaultLegendInterpreter','none');
    width=15.5;height=10;if k==1,height=8.5;elseif k==3,height=16;elseif k==4,height=11.5;end
    if k==1
        a=axes(fig,'Position',[.15,.17,.79,.75]);hold(a,'on');
        fine=linspace(0,max(259200,candidate.event.stop_time_s),1000)';
        gp=q4_compute("interpolate_radius",fine,in,"pchip");gl=q4_compute("interpolate_radius",fine,in,"linear");
        plot(a,fine/3600,gp.R_m*100,'Color',blue,'DisplayName','PCHIP');
        plot(a,fine/3600,gl.R_m*100,'--','Color',orange,'DisplayName','线性插值');
        plot(a,in.radius(:,1)/3600,in.radius(:,2),'.','Color',[.2,.2,.2],'MarkerSize',5,'DisplayName','给定半径');
        xline(a,72,':','72小时观测终点','HandleVisibility','off');
        xlabel(a,'时间（小时）');ylabel(a,'半径（厘米）');legend(a,'Location','northeast');
    elseif k==2
        a=axes(fig,'Position',[.13,.17,.68,.75]);hold(a,'on');
        rows=unique(round(linspace(1,numel(f.time_s),240)));xi=linspace(0,1,81)';
        C=zeros(numel(xi),numel(rows));R=f.radius_cm(rows)';times=f.time_s(rows)';
        for j=1:numel(rows)
            t=times(j);[found,index]=ismember(t,candidate.time_s);
            if found,y=candidate.states(:,index);
            else,assert(t==candidate.event.stop_time_s);y=candidate.event.stop_state;end
            g=q4_compute("interpolate_radius",t,in,p.geometry);bc=q4_compute("interpolate_ambient",t,in,p.environment);
            rec=reconstruct_q4_fields(y,candidate.mesh,g,bc,p,xi);C(:,j)=rec.C;
        end
        X=repmat(times/3600,numel(xi),1);Y=xi*R;
        surface(a,X,Y,zeros(size(C)),C,'EdgeColor','none','FaceColor','interp','HandleVisibility','off');
        view(a,2);colormap(a,parula(256));clim(a,[0,p.C0]);
        bar=colorbar(a,'Position',[.84,.17,.025,.75]);bar.Label.String='干基含水率（kg/kg）';
        plot(a,times/3600,R,'k-','LineWidth',1,'DisplayName','移动表面');
        contour(a,X,Y,C,[p.threshold,p.threshold],'LineColor',orange,'LineStyle','--', ...
            'LineWidth',1,'DisplayName','0.15等值线');
        xlabel(a,'时间（小时）');ylabel(a,'实际径向距离（厘米）');
        ylim(a,[0,p.R0_m*100]);xlim(a,[0,candidate.event.stop_time_h]);legend(a,'Location','northeast');
    elseif k==3
        a=axes(fig,'Position',[.15,.71,.79,.23]);hold(a,'on');
        plot(a,extra.progress(:,1)/3600,extra.progress(:,4),'Color',blue,'DisplayName','全域最大值');
        plot(a,extra.progress(:,1)/3600,extra.progress(:,2),'--','Color',orange,'DisplayName','截面平均');
        yline(a,.15,':','阈值0.15','HandleVisibility','off');ylabel(a,'含水率（kg/kg）');xlabel(a,'时间（小时）');legend(a,'Location','northeast');
        b=axes(fig,'Position',[.15,.40,.79,.21]);hold(b,'on');
        plot(b,extra.progress(:,1)/3600,extra.progress(:,3),'Color',green);
        ylim(b,[0,1.03]);xlabel(b,'时间（小时）');ylabel(b,'达标体积分数');
        if ~isempty(candidate.mean_event)
            xline(b,candidate.mean_event.root_s/3600,'--','平均首达','LabelVerticalAlignment','top');
        end
        xline(b,candidate.event.root_s/3600,':','全域临界','LabelVerticalAlignment','bottom');
        c=axes(fig,'Position',[.15,.11,.79,.17]);hold(c,'on');
        e=candidate.event;offset=[e.left_s,e.right_s]-e.root_s;
        plot(c,offset,[1,1],'-o','Color',blue,'DisplayName','候选根括区');
        plot(c,e.stop_time_s-e.root_s,1,'s','Color',orange,'DisplayName','候选停止点');
        xline(c,0,':','HandleVisibility','off');ylim(c,[.7,1.45]);
        set(c,'YTick',[]);xlabel(c,'相对候选临界时刻的秒数');
        text(c,.02,.86,'局部定位区间，不含网格与时间离散误差','Units','normalized','FontSize',7);
        legend(c,'Location','south','Orientation','horizontal','FontSize',6.5);
        text(a,0,1.06,'(a)','Units','normalized');text(b,0,1.06,'(b)','Units','normalized');
        text(c,0,1.06,'(c)','Units','normalized');
    else
        a=axes(fig,'Position',[.22,.63,.72,.28]);hold(a,'on');
        if isfield(summary,'cases')&&~isempty(summary.cases)
            records=summary.cases;
            labels=string({records.name});wanted=["main_N20","appendix4_fixed","appendix3_fixed","appendix3_shrinking"];
            readable=["附录4／收缩","附录4／固定","附录3／固定","附录3／收缩"];
            [found,ix]=ismember(wanted,labels);
            row=1;names=strings(1,0);
            for j=1:4
                if found(j)
                    plot(a,records(ix(j)).stop_h,row,'o','Color',blue);names(end+1)=readable(j);row=row+1; %#ok<AGROW>
                end
            end
            if row>1,set(a,'YTick',1:row-1,'YTickLabel',names);ylim(a,[.5,row-.5]);end
        end
        xlabel(a,'候选停止时间（小时）');ylabel(a,'物性/几何工况');
        b=axes(fig,'Position',[.15,.17,.79,.28]);hold(b,'on');
        rows=extra.mechanism(:,2)==0;
        z=extra.mechanism(rows,:);
        plot(b,z(:,1)/3600,z(:,8),'Color',blue,'DisplayName','含水率贡献');
        plot(b,z(:,1)/3600,z(:,9),'--','Color',orange,'DisplayName','温度贡献');
        plot(b,z(:,1)/3600,z(:,10),':','Color',green,'DisplayName','收缩贡献');
        xlabel(b,'时间（小时）');ylabel(b,'中心局部系数的对数变化');legend(b,'Location','best');
    end
    folder=fullfile(cfg.figuresDir,'q4',types(k));if ~isfolder(folder),mkdir(folder);end
    stem=fullfile(folder,stems(k));
    export_publication_figure(fig,stem,'Formats',"svg",'QADir',cfg.figureQADir, ...
        'Width',width,'Height',height,'Units',"centimeters",'DPI',300,'StrictDesign',true);
    q4_editable_svg(fig,stem+".svg",width,height);
    q4_normalize_svg_size(stem+".svg",width,height);
    paths(end+1)="figures/q4/"+types(k)+"/"+stems(k)+".svg";
end
end


function audit=q4_editable_svg(fig,destination,widthCm,heightCm)
fig.Units='centimeters';fig.Position(3:4)=[widthCm,heightCm];
fig.PaperUnits='centimeters';fig.PaperPosition=[0,0,widthCm,heightCm];
drawnow;texts=findall(fig,'Type','text','Visible','on');
keep=arrayfun(@(t)strlength(join(string(t.String)))>0,texts);texts=texts(keep);
fig.Units='pixels';sizePx=fig.Position(3:4);fig.Units='centimeters';
records=struct([]);
for k=1:numel(texts)
    t=texts(k);assert(string(t.Interpreter)=="none",'CUMCM:q4:SvgText','Editable adapter requires plain text.');
    old=t.Units;t.Units='pixels';pos=t.Position;t.Units=old;
    ax=ancestor(t,'axes');assert(~isempty(ax),'CUMCM:q4:SvgText','Unsupported text parent.');
    box=getpixelposition(ax,true);record=struct('x',(box(1)+pos(1)-1)/sizePx(1), ...
        'y',1-(box(2)+pos(2)-1)/sizePx(2),'text',string(t.String), ...
        'font',string(t.FontName),'size',t.FontSize,'weight',string(t.FontWeight), ...
        'color',t.Color,'rotation',t.Rotation,'align',string(t.HorizontalAlignment), ...
        'valign',string(t.VerticalAlignment));
    if isempty(records),records=record;else,records(k)=record;end
    t.Visible='off';
end
restore=onCleanup(@()set(texts,'Visible','on'));
native=string(tempname(fileparts(destination)))+".svg";print(fig,char(native),'-dsvg','-painters','-r600');
xml=fileread(native);
pixelSizes=[];
w=regexp(xml,'<svg[^>]*\swidth="([\d.]+)"','tokens','once');
h=regexp(xml,'<svg[^>]*\sheight="([\d.]+)"','tokens','once');
assert(~isempty(w)&&~isempty(h),'CUMCM:q4:SvgText','SVG coordinate dimensions missing.');
width=str2double(w{1});height=str2double(h{1});pieces=cell(numel(records),1);
for k=1:numel(records)
    r=records(k);anchor='start';if r.align=="center",anchor='middle';elseif r.align=="right",anchor='end';end
    dy=0;if r.valign=="top"||r.valign=="cap",dy=0.8;elseif r.valign=="middle",dy=0.35;elseif r.valign=="bottom",dy=-0.15;end
    x=width*r.x;y=height*r.y;fontSize=r.size*width/(widthCm/2.54*72);
    pieces{k}=sprintf(['<text x="%.8g" y="%.8g" dy="%.3gem" text-anchor="%s" ', ...
        'transform="rotate(%.8g %.8g %.8g)" style="font-family:%s;font-size:%.8gpx;', ...
        'font-weight:%s;fill:rgb(%d,%d,%d);stroke:none">%s</text>'], ...
        x,y,dy,anchor,-r.rotation,x,y,escape(r.font),fontSize,r.weight,round(255*r.color),escape(r.text));
end
xml=regexprep(xml,'</svg\s*>',['<g id="editable-matlab-labels">',strjoin(pieces,''),'</g></svg>']);
xml=regexprep(xml,'(<svg[^>]*\swidth=")[\d.]+(")',sprintf('$1%.8gcm$2',widthCm),'once');
xml=regexprep(xml,'(<svg[^>]*\sheight=")[\d.]+(")',sprintf('$1%.8gcm$2',heightCm),'once');
xml=regexprep(xml,'<svg ',sprintf('<svg viewBox="0 0 %.8g %.8g" ',width,height),'once');
fid=fopen(destination,'w','n','UTF-8');assert(fid>=0);clean=onCleanup(@()fclose(fid));fprintf(fid,'%s',xml);clear clean
qaFolder=fullfile(fileparts(fileparts(fileparts(fileparts(destination)))),'artifacts','q4','svg_native');
if ~isfolder(qaFolder),mkdir(qaFolder);end
[~,stem]=fileparts(destination);movefile(native,fullfile(qaFolder,string(stem)+"_no_text.svg"),'f');
previewDir=fullfile(fileparts(fileparts(qaFolder)),'validation','figure_previews');
if ~isfolder(previewDir),mkdir(previewDir);end
transcoder=javaObject('org.apache.batik.transcoder.image.PNGTranscoder');
transcoder.addTranscodingHint(org.apache.batik.transcoder.SVGAbstractTranscoder.KEY_WIDTH,java.lang.Float(widthCm/2.54*300));
transcoder.addTranscodingHint(org.apache.batik.transcoder.SVGAbstractTranscoder.KEY_HEIGHT,java.lang.Float(heightCm/2.54*300));
cleanXml=regexprep(xml,'<!DOCTYPE[\s\S]*?>','');
input=javaObject('org.apache.batik.transcoder.TranscoderInput',java.io.StringReader(cleanXml));
renderPath=fullfile(previewDir,string(stem)+"_color.png");stream=java.io.FileOutputStream(char(renderPath));
streamCleanup=onCleanup(@()stream.close());
transcoder.transcode(input,javaObject('org.apache.batik.transcoder.TranscoderOutput',stream));clear streamCleanup
rgb=imread(renderPath);ppm=round(300/.0254);
imwrite(rgb,renderPath,'png','ResolutionUnit','meter','XResolution',ppm,'YResolution',ppm);
gray=uint8(round(.298936*double(rgb(:,:,1))+.587043*double(rgb(:,:,2))+.114021*double(rgb(:,:,3))));
imwrite(gray,fullfile(previewDir,string(stem)+"_gray.png"),'png','ResolutionUnit','meter','XResolution',ppm,'YResolution',ppm);
audit=struct('visible_text_count',numel(texts),'serialized_text_count',numel(records), ...
    'width_cm',widthCm,'height_cm',heightCm,'native_matlab_data_geometry_preserved',true,'native_raster_export_dpi',600, ...
    'embedded_image_pixels_width_height',pixelSizes,'preview_from_final_svg',true);
end
function s=escape(s)
s=replace(string(s),["&","<",">",'"'],["&amp;","&lt;","&gt;","&quot;"]);
end


function q4_normalize_svg_size(file,width_cm,height_cm)
text=fileread(file);header=regexp(text,'<svg\s[^>]*>','match','once');assert(~isempty(header));
updated=header;
if ~contains(updated,'viewBox=')
    canvas=regexp(text,'<path d="M0 0 L([0-9.]+) 0 L[0-9.]+ ([0-9.]+) L0 ','tokens','once');
    assert(numel(canvas)==2,'CUMCM:q4:SvgCanvas','Cannot identify original SVG canvas.');
    updated=regexprep(updated,'>$',sprintf(' viewBox="0 0 %s %s">',canvas{1},canvas{2}));
end
updated=regexprep(updated,'\swidth="[^"]*"',sprintf(' width="%.15gcm"',width_cm),'once');
updated=regexprep(updated,'\sheight="[^"]*"',sprintf(' height="%.15gcm"',height_cm),'once');
assert(~strcmp(header,updated)||contains(header,sprintf('width="%.15gcm"',width_cm)));
text=replace(text,header,updated);fid=fopen(file,'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',text); %#ok<NASGU>
end


function q4_pack_ooxml(parts,destination)
listing=dir(fullfile(parts,'**','*')); listing=listing(~[listing.isdir]);
relative=strings(numel(listing),1);
for k=1:numel(listing)
    relative(k)=replace(extractAfter(string(fullfile(listing(k).folder,listing(k).name)),strlength(string(parts))+1),filesep,'/');
end
[relative,order]=sort(relative); listing=listing(order);
out=java.util.zip.ZipOutputStream(java.io.FileOutputStream(char(destination)));
clean=onCleanup(@()out.close()); out.setLevel(6);
for k=1:numel(listing)
    assert(~startsWith(relative(k),'/') && ~contains(relative(k),'../'));
    entry=java.util.zip.ZipEntry(char(relative(k))); entry.setTime(946684800000);
    out.putNextEntry(entry);
    bytes=java.nio.file.Files.readAllBytes(java.io.File(fullfile(listing(k).folder,listing(k).name)).toPath());
    out.write(bytes,0,numel(bytes)); out.closeEntry();
end
end


function outputDir=export_q4_boundary_deliverables(cfg)
projectRoot=cfg.projectRoot;p=q4_compute("parameters",cfg);in=q4_compute("load_inputs",cfg);
stats=in.environment_selection;q4_validation("environment",in,p);
outputDir=fullfile(projectRoot,'deliverables','q4_统一边界代码图表');
if ~isfolder(outputDir),mkdir(outputDir);end
figureDir=fullfile(outputDir,'图片');if ~isfolder(figureDir),mkdir(figureDir);end
q={stats.temperature,stats.moisture};labels=["温度","水分浓度"];units=["℃","kg/kg"];
summary=zeros(2,12);
book=fullfile(outputDir,'边界选窗结果.xlsx');
for c=1:2
    s=q{c}.selection;w=q{c}.windows;
    summary(c,:)=[s.first_exceed_min,s.start_min,s.end_min,s.duration_min,s.sample_count, ...
        s.mean_value,s.sample_variance,s.variance_penalty,s.sample_penalty,s.total_loss, ...
        s.candidate_count,s.screen_pass_count];
    values=[w.start_min,w.end_min,w.duration_min,w.sample_count,w.mean_value,w.sample_variance, ...
        w.block_sizes,w.block_means,w.block_variances,w.standardized_drift,w.three_block_mean_span, ...
        w.three_block_variance_ratio,w.variance_penalty,w.sample_penalty,w.total_loss];
    headers=["窗口起点_min","窗口终点_min","时长_min","样本数","算术均值","样本方差", ...
        "前段样本数","中段样本数","后段样本数","前段均值","中段均值","后段均值", ...
        "前段方差","中段方差","后段方差","漂移指标D","均值跨度B","方差比R", ...
        "加权方差惩罚","加权样本惩罚","总损失L"];
    tab=array2table(values,'VariableNames',headers);
    tab.("通过稳定筛选")=w.screen_pass;tab.("选中")=w.selected;tab.("未通过原因")=w.rejection_reason;
    writetable(tab,fullfile(outputDir,labels(c)+'候选窗口.csv'),'Encoding','UTF-8');
    writetable(tab,book,'Sheet',labels(c)+'候选','WriteMode','overwritesheet');
    drawFigure(c,w,s);
end
tab=array2table(summary,'VariableNames',["首次越阈_min","窗口起点_min","窗口终点_min", ...
    "时长_min","样本数","平衡水平估计","样本方差","加权方差惩罚","加权样本惩罚", ...
    "总损失L","候选数","通过筛选数"]);
tab=addvars(tab,labels',units','Before',1,'NewVariableNames',{'信号','均值单位'});
writetable(tab,fullfile(outputDir,'边界选值汇总.csv'),'Encoding','UTF-8');
writetable(tab,book,'Sheet','边界选值汇总','WriteMode','overwritesheet');
notes={"计算口径","首次严格越阈后的全部后缀；至少30分钟；均值为算术平均；方差分母n-1"; ...
    "损失函数","0.8 ln(方差/全部候选最小方差)+0.2 ln(最多样本数/窗口样本数)"; ...
    "筛选条件","D≤0.5，B≤0.5，R≤3；不合格窗口不参与最小化"; ...
    "方差单位","温度为℃平方；水分浓度为(kg/kg)平方"; ...
    "边界使用","0至4小时实测插值不改；仅4小时以后同时使用两组选中均值"; ...
    "解释范围","窗口起点不等于真实物理平衡时刻；本表不包含新干燥时间"};
writecell(notes,book,'Sheet','说明','WriteMode','overwritesheet');
write_q4_environment_evidence(cfg,in,'deliverables/q4_统一边界代码图表');
for c=1:2
    back=readmatrix(book,'Sheet',labels(c)+'候选','Range','A2:U'+string(height(q{c}.windows)+1),'UseExcel',false);
    assert(size(back,1)==height(q{c}.windows)&&all(isfinite(back),'all'));
end
fprintf('Boundary deliverables saved: %s\n',outputDir);

    function drawFigure(c,w,s)
        fig=figure('Visible','off','Color','w');cleanup=onCleanup(@()close(fig));
        apply_publication_style(fig,'en','report');
        set(fig,'Units','centimeters','Position',[2,2,18,16], ...
            'DefaultAxesFontName','Microsoft YaHei UI','DefaultTextFontName','Microsoft YaHei UI', ...
            'DefaultAxesFontSize',9,'DefaultTextInterpreter','none','DefaultLegendInterpreter','none');
        blue=[0,.447,.698];orange=[.835,.369,0];gray=[.52,.52,.52];green=[0,.45,.32];
        ax=axes(fig,'Position',[.15,.59,.81,.31]);hold(ax,'on');
        first=s.first_exceed_min;t=in.ambient(:,1)/60;mask=t>=first;y=in.ambient(mask,c+1);
        span=max(y)-min(y);yl=[min(y)-.12*span,max(y)+.12*span];
        patch(ax,[s.start_min,240,240,s.start_min],[yl(1),yl(1),yl(2),yl(2)], ...
            [.90,.95,.92],'EdgeColor','none','HandleVisibility','off');
        plot(ax,t(mask),y,'Color',blue,'LineWidth',.9,'DisplayName','原始观测');
        plot(ax,[s.start_min,240],[s.mean_value,s.mean_value],'--','Color',orange, ...
            'LineWidth',1.6,'DisplayName','选中窗口均值');
        xline(ax,s.start_min,':','Color',green,'LineWidth',1.3,'DisplayName','窗口起点');
        yline(ax,s.threshold,':','Color',gray,'DisplayName','搜索阈值');
        xlim(ax,[first-1,241]);ylim(ax,yl);xlabel(ax,'观测时间（分钟）');ylabel(ax,labels(c)+'（'+units(c)+'）');
        if c==2,ytickformat(ax,'%.5f');end
        legend(ax,'Location','northoutside','Orientation','horizontal','FontSize',8,'NumColumns',4);
        text(ax,0,1.29,sprintf('(a) 选中%d—240分钟，共%d点；阴影为取样区间',s.start_min,s.sample_count), ...
            'Units','normalized','FontSize',9);
        ax2=axes(fig,'Position',[.15,.13,.81,.27]);hold(ax2,'on');
        plot(ax2,w.start_min,w.variance_penalty,'--','Color',blue,'DisplayName','方差惩罚');
        plot(ax2,w.start_min,w.sample_penalty,':','Color',orange,'LineWidth',1.6,'DisplayName','样本量惩罚');
        plot(ax2,w.start_min,w.total_loss,'-','Color',[.15,.15,.15],'DisplayName','总损失');
        plot(ax2,w.start_min(~w.screen_pass),w.total_loss(~w.screen_pass),'x','Color',gray, ...
            'LineStyle','none','MarkerSize',5,'DisplayName','未通过筛选');
        plot(ax2,s.start_min,s.total_loss,'o','Color',green,'MarkerFaceColor',green, ...
            'MarkerSize',6,'DisplayName','选中窗口');
        xlabel(ax2,'候选窗口起点（分钟）');ylabel(ax2,'加权惩罚 / 总损失（无量纲）');
        xlim(ax2,[first-1,211]);ylim(ax2,[0,max(w.total_loss)*1.12]);
        legend(ax2,'Location','northoutside','Orientation','horizontal','NumColumns',3,'FontSize',8);
        text(ax2,0,1.38,'(b) 仅在通过稳定筛选的窗口中最小化总损失', ...
            'Units','normalized','FontSize',9);
        stem=fullfile(figureDir,labels(c)+'选窗');
        export_publication_figure(fig,stem,'Formats',"svg",'QADir',figureDir, ...
            'Width',18,'Height',16,'Units',"centimeters",'DPI',300,'StrictDesign',true);
        q4_editable_svg(fig,stem+'.svg',18,16);q4_normalize_svg_size(stem+'.svg',18,16);
        savefig(fig,stem+'.fig');
    end
end


function style = apply_publication_style(fig, language, widthProfile)
arguments
    fig (1,1) matlab.ui.Figure
    language (1,1) string {mustBeMember(language,["zh","en"])} = "zh"
    widthProfile (1,1) string {mustBeMember(widthProfile,["single","double","report"])} = "report"
end

widths = struct("single", 3.5, "double", 7.2, "report", 6.3);
fontName = chooseFont(language);
colors = [
    0.0000 0.4471 0.6980
    0.9020 0.6235 0.0000
    0.0000 0.6196 0.4510
    0.8353 0.3686 0.0000
    0.8000 0.4745 0.6549
    0.3373 0.7059 0.9137
    0.4196 0.4471 0.5020
];

widthIn = widths.(widthProfile);
fig.Units = "inches";
fig.Position(3:4) = [widthIn, widthIn * 0.62];
fig.Color = "white";
set(fig, "DefaultAxesFontName", fontName, ...
    "DefaultAxesFontSize", 7.5, ...
    "DefaultAxesLineWidth", 0.7, ...
    "DefaultAxesTitleFontSizeMultiplier", 1.0, ...
    "DefaultAxesTitleFontWeight", "normal", ...
    "DefaultAxesLabelFontSizeMultiplier", 1.0, ...
    "DefaultAxesColorOrder", colors, ...
    "DefaultLineLineWidth", 1.1, ...
    "DefaultLineMarkerSize", 3.5, ...
    "DefaultLegendBox", "off");

axesList = findall(fig, "Type", "axes");
for k = 1:numel(axesList)
    axesList(k).FontName = fontName;
    axesList(k).FontSize = 7.5;
    axesList(k).LineWidth = 0.7;
    axesList(k).TitleFontSizeMultiplier = 1.0;
    axesList(k).TitleFontWeight = "normal";
    axesList(k).LabelFontSizeMultiplier = 1.0;
    axesList(k).ColorOrder = colors;
    axesList(k).Box = "off";
    grid(axesList(k), "off");
end

style = struct("font", fontName, "colors", colors, ...
    "sizeInches", [widthIn, widthIn * 0.62]);
end

function fontName = chooseFont(language)
fonts = string(listfonts);
if language == "zh"
    candidates = ["Noto Sans CJK SC","Source Han Sans SC", ...
        "Microsoft YaHei","SimHei","PingFang SC"];
else
    candidates = ["Arial","Helvetica","Times New Roman"];
end
fontName = "Helvetica";
for candidate = candidates
    if any(strcmpi(fonts, candidate))
        fontName = candidate;
        return;
    end
end
warning("未找到首选字体，导出后必须检查中文与特殊符号是否缺字。");
end


function issues = audit_publication_figure(fig)
arguments
    fig (1,1) matlab.ui.Figure
end

drawnow;
issues = strings(0,1);
axesList = findall(fig, "Type", "axes");

for k = 1:numel(axesList)
    ax = axesList(k);
    originalUnits = ax.Units;
    ax.Units = "normalized";
    position = ax.Position;
    inset = ax.TightInset;
    contentBounds = [
        position(1) - inset(1)
        position(2) - inset(2)
        position(1) + position(3) + inset(3)
        position(2) + position(4) + inset(4)
    ];
    ax.Units = originalUnits;
    if any(contentBounds(1:2) < -0.005) || any(contentBounds(3:4) > 1.005)
        issues(end+1) = "第 " + k + ...
            " 个坐标轴的标签或刻度可能超出画布"; %#ok<AGROW>
    end

    titleParts = string(ax.Title.String);
    titleText = strtrim(join(titleParts(:), " "));
    if strlength(titleText) > 30
        issues(end+1) = "第 " + k + ...
            " 个坐标轴标题超过 30 字符；只允许必要的短面板标题，完整说明写入论文题注"; %#ok<AGROW>
    end

    lineList = findall(ax, "Type", "line");
    for line = reshape(lineList, 1, [])
        marker = string(line.Marker);
        if marker == "none" || marker == ""
            continue;
        end
        hasSparseMarkers = isprop(line, "MarkerIndices") && ~isempty(line.MarkerIndices);
        if numel(line.XData) > 25 && ~hasSparseMarkers
            issues(end+1) = "第 " + k + ...
                " 个坐标轴对稠密曲线逐点绘制标记；应取消标记或设置 MarkerIndices"; %#ok<AGROW>
            break;
        end
    end

    barList = findall(ax, "Type", "bar");
    for barObject = reshape(barList, 1, [])
        isHorizontal = isprop(barObject, "Horizontal") && ...
            string(barObject.Horizontal) == "on";
        if isHorizontal
            limits = ax.XLim;
        else
            limits = ax.YLim;
        end
        tolerance = max(abs(diff(limits)), 1) * 1e-9;
        if abs(barObject.BaseValue) > tolerance || ...
                limits(1) > tolerance || limits(2) < -tolerance
            issues(end+1) = "第 " + k + ...
                " 个坐标轴的柱状图未从零开始；若必须截断，应改用点图/区间图并显式说明"; %#ok<AGROW>
            break;
        end
    end

    imageList = findall(ax, "Type", "image");
    textList = findall(ax, "Type", "text");
    hasOwnColorbar = false;
    if isprop(ax, "Colorbar") && ~isempty(ax.Colorbar)
        hasOwnColorbar = isvalid(ax.Colorbar);
    elseif isappdata(ax, "ColorbarPeerHandle")
        colorbarHandle = getappdata(ax, "ColorbarPeerHandle");
        hasOwnColorbar = ~isempty(colorbarHandle) && isvalid(colorbarHandle);
    end
    for imageObject = reshape(imageList, 1, [])
        valueCount = numel(imageObject.CData);
        if valueCount <= 4 && numel(textList) >= valueCount && hasOwnColorbar
            issues(end+1) = "第 " + k + ...
                " 个坐标轴是已标数值的 2×2 小矩阵，不应再使用冗余 colorbar"; %#ok<AGROW>
            break;
        end
    end
end

legendList = findall(fig, "Type", "legend");
for legendObject = reshape(legendList, 1, [])
    labels = string(legendObject.String);
    if numel(labels) > 5
        issues(end+1) = "图例超过 5 项；应直接标注、改为共享图例或拆分证据"; %#ok<AGROW>
        break;
    end
end

issues = unique(issues, "stable");
end


function outputs = export_publication_figure(fig, outputStem, options)
arguments
    fig (1,1) matlab.ui.Figure
    outputStem (1,1) string
    options.Formats (1,:) string = "svg"
    options.DPI (1,1) double {mustBeInteger,mustBePositive} = 300
    options.QADir (1,1) string = ""
    options.StrictDesign (1,1) logical = true
    options.Width (1,1) double {mustBePositive} = 6.3
    options.Height (1,1) double {mustBePositive} = 3.9
    options.Units (1,1) string {mustBeMember(options.Units, ...
        ["inches","centimeters","points"])} = "inches"
    options.Padding (1,1) string {mustBeMember(options.Padding, ...
        ["tight","figure"])} = "figure"
end

formats = unique(lower(options.Formats), "stable");
if ~isequal(formats, "svg")
    error("CUMCM:UnsupportedFigureFormat", ...
        "正式交付接口只允许 Formats=""svg""；PNG 仅由 QADir 生成内部预览。");
end
[folder, name, ~] = fileparts(outputStem);
if strlength(folder) == 0
    folder = ".";
end
if ~isfolder(folder)
    mkdir(folder);
end

drawnow;
designIssues = audit_publication_figure(fig);
if options.StrictDesign && ~isempty(designIssues)
    error("CUMCM:FigureDesignAuditFailed", ...
        "出版设计预检未通过：%s", join(designIssues, "；"));
end

originalUnits = fig.Units;
originalPosition = fig.Position;
originalPaperUnits = fig.PaperUnits;
originalPaperPosition = fig.PaperPosition;
restore = onCleanup(@() restoreFigure(fig, originalUnits, originalPosition, ...
    originalPaperUnits, originalPaperPosition)); %#ok<NASGU>
fig.Units = options.Units;
position = fig.Position;
position(3:4) = [options.Width, options.Height];
fig.Position = position;
fig.PaperUnits = options.Units;
fig.PaperPosition = [0, 0, options.Width, options.Height];
drawnow;

supportsSizedExportgraphics = ~verLessThan("matlab", "25.1");
common = {"Width", options.Width, "Height", options.Height, ...
    "Units", options.Units, "Padding", options.Padding};
formalPaths = strings(0,1);
for format = formats
    pathValue = fullfile(folder, name + "." + format);
    if supportsSizedExportgraphics
        exportgraphics(fig, pathValue, common{:}, "ContentType", "vector");
    else
        print(fig, char(pathValue), "-dsvg", "-painters");
    end
    formalPaths(end+1,1) = pathValue; %#ok<AGROW>
end

qaColor = "";
qaGray = "";
if strlength(options.QADir) > 0
    if ~isfolder(options.QADir)
        mkdir(options.QADir);
    end
    qaColor = fullfile(options.QADir, name + "_color.png");
    qaGray = fullfile(options.QADir, name + "_gray.png");
    if supportsSizedExportgraphics
        exportgraphics(fig, qaColor, common{:}, "Resolution", options.DPI);
        exportgraphics(fig, qaGray, common{:}, "Resolution", options.DPI, ...
            "Colorspace", "gray");
    else
        print(fig, char(qaColor), "-dpng", char("-r" + options.DPI));
        writeGrayscalePreview(qaColor, qaGray, options.DPI);
    end
end

outputs = struct("formal", formalPaths, ...
    "qa_color", qaColor, "qa_grayscale", qaGray, ...
    "dpi", options.DPI, "width", options.Width, ...
    "height", options.Height, "units", options.Units, ...
    "backend", conditionalBackend(supportsSizedExportgraphics));
end

function value = conditionalBackend(supported)
if supported
    value = "exportgraphics-sized";
else
    value = "print-compatibility";
end
end

function writeGrayscalePreview(colorPath, grayPath, dpi)
rgb = imread(colorPath);
if ndims(rgb) == 3 && size(rgb, 3) >= 3
    gray = uint8(round(0.298936 * double(rgb(:,:,1)) + ...
        0.587043 * double(rgb(:,:,2)) + 0.114021 * double(rgb(:,:,3))));
else
    gray = rgb;
end
pixelsPerMeter = round(dpi / 0.0254);
imwrite(gray, char(grayPath), "png", "ResolutionUnit", "meter", ...
    "XResolution", pixelsPerMeter, "YResolution", pixelsPerMeter);
end

function restoreFigure(fig, units, position, paperUnits, paperPosition)
if isvalid(fig)
    fig.Units = units;
    fig.Position = position;
    fig.PaperUnits = paperUnits;
    fig.PaperPosition = paperPosition;
end
end
