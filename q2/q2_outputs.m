function varargout=q2_outputs(action,varargin)
action=string(action);
assert(isscalar(action),'CUMCM:q2:OutputAction','Invalid output action.');
switch action
    case "validation_io"
        value=q2_validation_io(varargin{:});
        if nargout>0, varargout{1}=value; end
    case "write_csv"
        varargout{1}=write_q2_csv(varargin{:});
    case "write_result2"
        varargout{1}=write_q2_result2(varargin{:});
    case "write_restart"
        varargout{1}=write_q2_restart(varargin{:});
    case "plot_results"
        varargout{1}=plot_q2_results(varargin{:});
    otherwise
        error('CUMCM:q2:OutputAction','Unknown output action: %s',action);
end
end

function value=q2_validation_io(action,path,data)
switch string(action)
    case "json"
        folder=fileparts(path);if ~isfolder(folder),mkdir(folder);end
        fid=fopen(path,'w','n','UTF-8');assert(fid>=0,'CUMCM:q2:OutputIO','Cannot write JSON.');
        clean=onCleanup(@()fclose(fid));
        fprintf(fid,'%s\n',jsonencode(data,'PrettyPrint',true));value=path;
    case "hash"
        md=java.security.MessageDigest.getInstance('SHA-256');
        stream=java.io.FileInputStream(char(path));clean=onCleanup(@()stream.close());
        channel=stream.getChannel();buffer=java.nio.ByteBuffer.allocate(1024^2);
        while channel.read(buffer)>0
            buffer.flip();md.update(buffer);buffer.clear();
        end
        value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
    case "record"
        info=dir(fullfile(path,data));assert(isscalar(info),'CUMCM:q2:MissingEvidence','Missing artifact.');
        value=struct('path',string(data),'sha256',q2_validation_io('hash',fullfile(path,data)), ...
            'bytes',info.bytes);
    otherwise
        error('CUMCM:q2:OutputIO','Unknown evidence operation.');
end
end

function paths=write_q2_csv(cfg,fields,mechanism,validation)
folder=fullfile(cfg.resultsDir,'q2');if ~isfolder(folder),mkdir(folder);end
paths="results/q2/"+["q2_temperature_unrounded.csv","q2_moisture_unrounded.csv", ...
    "q2_table3_temperature.csv","q2_table4_moisture.csv","q2_mechanism.csv","q2_validation_summary.json"];
names={'temperature_C','moisture_kgkg'};
for field=1:2
    data=fields.(names{field});assert(isequal(size(data),[10800,21])&&all(isfinite(data),'all'));
    fid=fopen(fullfile(cfg.projectRoot,paths(field)),'w','n','UTF-8');assert(fid>=0);
    clean=onCleanup(@()fclose(fid));fprintf(fid,'时间\\到药材中心的距离');
    fprintf(fid,',%.1f',fields.radius_cm);fprintf(fid,'\n');
    fprintf(fid,['%.0f',repmat(',%.17g',1,21),'\n'],[fields.time_s,data]');clear clean
    readback=readmatrix(fullfile(cfg.projectRoot,paths(field)),'NumHeaderLines',1);
    assert(isequal(readback,[fields.time_s,data]),'CUMCM:q2:OutputMapping','Full CSV roundtrip failed.');
    fid=fopen(fullfile(cfg.projectRoot,paths(field+2)),'w','n','UTF-8');assert(fid>=0);
    clean=onCleanup(@()fclose(fid));fprintf(fid,'时间（小时）');
    fprintf(fid,',半径%.1f厘米',fields.radius_cm(fields.paper_columns));fprintf(fid,'\n');
    fprintf(fid,['%.1f',repmat(',%.4f',1,5),'\n'], ...
        [fields.paper_time_s/3600,round(data(fields.paper_time_s,fields.paper_columns),4)]');clear clean
    paper=readmatrix(fullfile(cfg.projectRoot,paths(field+2)),'NumHeaderLines',1);
    assert(isequal(paper(:,2:end),round(data(fields.paper_time_s,fields.paper_columns),4)));
end
fid=fopen(fullfile(cfg.projectRoot,paths(5)),'w','n','UTF-8');assert(fid>=0);
clean=onCleanup(@()fclose(fid));
fprintf(fid,'时间（秒）,位置,温度（℃）,干基含水率（kg/kg）,密度（kg/m³）,比热容（J/(kg·K)）,导热系数（W/(m·K)）,水分扩散系数（m²/s）,热扩散率（m²/s）,含水率对对数扩散变化的贡献（无量纲）,温度对对数扩散变化的贡献（无量纲）\n');
for position=1:2
    label=["中心","表面"];
    numeric=[mechanism.temperature_C(:,position),mechanism.moisture_kgkg(:,position), ...
        mechanism.rho(:,position),mechanism.cp(:,position),mechanism.k(:,position), ...
        mechanism.D(:,position),mechanism.alpha(:,position),mechanism.logD_C(:,position),mechanism.logD_T(:,position)];
    for row=1:10800
        fprintf(fid,'%.0f,%s',fields.time_s(row),label(position));fprintf(fid,',%.17g',numeric(row,:));fprintf(fid,'\n');
    end
end
clear clean
q2_validation_io('json',fullfile(cfg.projectRoot,paths(6)),validation);
end

function path=write_q2_result2(cfg,fields)
source=fullfile(cfg.rawDataDir,'附件3','result2.xlsx'); originalHash=hashFile(source);
folder=fullfile(cfg.processedDataDir,'q2'); if ~isfolder(folder), mkdir(folder); end
work=fullfile(cfg.projectRoot,'artifacts','q2','output_work'); if ~isfolder(work), mkdir(work); end
parts=string(tempname(work)); mkdir(parts);
copy=string(tempname(work))+".xlsx"; copyfile(source,copy); unzip(copy,parts);
stylesPath=fullfile(parts,'xl','styles.xml'); styles=fileread(stylesPath);
token=regexp(styles,'<cellXfs count="4">(.*?)</cellXfs>','tokens','once'); assert(~isempty(token));
xfs=regexp(token{1},'<xf\s[^>]*>.*?</xf>','match'); assert(numel(xfs)==4);
temperatureStyle=replace(xfs{3},'numFmtId="0"','numFmtId="164"');
moistureStyle=replace(xfs{2},'numFmtId="0"','numFmtId="164"');
temperatureStyle=replace(temperatureStyle,'<xf ','<xf applyNumberFormat="1" ');
moistureStyle=replace(moistureStyle,'<xf ','<xf applyNumberFormat="1" ');
styles=regexprep(styles,'(<styleSheet[^>]*>)', ...
    '$1<numFmts count="1"><numFmt numFmtId="164" formatCode="0.0000"/></numFmts>','once');
styles=replace(styles,['<cellXfs count="4">',token{1},'</cellXfs>'], ...
    ['<cellXfs count="6">',token{1},temperatureStyle,moistureStyle,'</cellXfs>']);
writeText(stylesPath,styles);
values={round(fields.temperature_C,4),round(fields.moisture_kgkg,4)};
for sheet=1:2
    pathXml=fullfile(parts,'xl','worksheets',"sheet"+sheet+".xml"); xml=fileread(pathXml);
    assert(~contains(xml,'<mergeCells')&&~contains(xml,'<f>'));
    rows=cell(10801,1); header=cell(1,22);
    header{1}='<c r="A1" s="2" t="s"><v>1</v></c>';
    headStyle=3-sheet; dataStyle=3+sheet;
    for col=2:22
        header{col}=sprintf('<c r="%s1" s="%d"><v>%.1f</v></c>',char(64+col),headStyle,fields.radius_cm(col-1));
    end
    rows{1}=['<row r="1" spans="1:22" ht="14.1" customHeight="1">',strjoin(header,''),'</row>'];
    for time=1:10800
        cells=cell(1,22); row=time+1;
        cells{1}=sprintf('<c r="A%d" s="1"><v>%d</v></c>',row,time);
        for col=2:22
            cells{col}=sprintf('<c r="%s%d" s="%d"><v>%.4f</v></c>',char(64+col),row,dataStyle,values{sheet}(time,col-1));
        end
        rows{row}=sprintf('<row r="%d" spans="1:22" ht="14.1" customHeight="1">%s</row>',row,strjoin(cells,''));
    end
    xml=regexprep(xml,'<sheetData>.*?</sheetData>',['<sheetData>',strjoin(rows,''),'</sheetData>']);
    xml=replace(xml,'ref="A1:F5"','ref="A1:V10801"');
    xml=replace(xml,'min="2" max="6"','min="2" max="22"');
    xml=replace(xml,'min="7" max="16384"','min="23" max="16384"');
    writeText(pathXml,xml);
end
shared=fullfile(parts,'xl','sharedStrings.xml'); ss=fileread(shared);
ss=regexprep(ss,' count="\d+"',' count="2"','once'); writeText(shared,ss);
first=string(tempname(work))+".xlsx"; second=string(tempname(work))+".xlsx";
q2_pack_ooxml(parts,first); canonicalize_q2_xlsx(first,cfg.projectRoot);
q2_pack_ooxml(parts,second); canonicalize_q2_xlsx(second,cfg.projectRoot);
assert(strcmp(hashFile(first),hashFile(second)),'CUMCM:q2:NonDeterministicArtifact','Workbook writes differ.');
sheets=sheetnames(first); assert(isequal(sheets(:),["温度";"水分浓度"]));
for k=1:2
    top=readcell(first,'Sheet',sheets(k),'Range','A1','UseExcel',false);
    assert(string(top{1})=="时间\到药材中心的距离");
    matrix=readmatrix(first,'Sheet',sheets(k),'Range','A1:V10801','UseExcel',false);
    assert(all(abs(matrix(1,2:end)-fields.radius_cm)<1e-12) && isequal(matrix(2:end,1),fields.time_s));
    assert(isequal(matrix(2:end,2:end),values{k}),'CUMCM:q2:OutputMappingFailure','XLSX values mismatch.');
end
assert(strcmp(hashFile(source),originalHash),'CUMCM:q2:TemplateMutation','Raw template changed.');
output=fullfile(folder,'result2.xlsx');
if isfile(output), copyfile(output,string(tempname(work))+"_previous.xlsx"); end
movefile(first,output,'f'); path="data/processed/q2/result2.xlsx";
end
function writeText(path,text)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',text);
end
function value=hashFile(path)
bytes=java.nio.file.Files.readAllBytes(java.io.File(char(path)).toPath());
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end

function canonicalize_q2_xlsx(path,projectRoot)
root=string(java.io.File(char(projectRoot)).getCanonicalPath());
actual=string(java.io.File(char(path)).getCanonicalPath());
assert(startsWith(lower(actual),lower(root+filesep)) && ~contains(lower(actual),lower(filesep+"data"+filesep+"raw"+filesep)), ...
    'CUMCM:q2:TemplateMutation','Only project-generated workbooks may be canonicalized.');
work=fullfile(root,'artifacts','q2','output_work'); if ~isfolder(work), mkdir(work); end
parts=string(tempname(work)); mkdir(parts); unzip(path,parts);
core=fullfile(parts,'docProps','core.xml');
if isfile(core)
    xml=fileread(core);
    for tag={'created','modified'}
        pattern=['(<dcterms:',tag{1},'[^>]*>)([^<]*)(</dcterms:',tag{1},'>)'];
        token=regexp(xml,pattern,'tokens','once');
        if ~isempty(token)
            xml=replace(xml,strjoin(token,''),[token{1},'2000-01-01T00:00:00Z',token{3}]);
        end
    end
    fid=fopen(core,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',xml); clear c
end
q2_pack_ooxml(parts,path);
end

function q2_pack_ooxml(parts,destination)
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

function paths=write_q2_restart(cfg,solution,bindings)
folder=fullfile(cfg.resultsDir,'q2');if ~isfolder(folder),mkdir(folder);end
assert(solution.final_time_s==10800&&numel(solution.final_state)==2*numel(solution.radius_m));
final_state=struct('time_s',solution.final_time_s,'radius_m',solution.radius_m, ...
    'temperature_K',solution.final_state(1:2:end),'moisture_kgkg',solution.final_state(2:2:end));
paths=["results/q2/q2_restart_10800s.mat","results/q2/q2_restart_metadata.json"];
save(fullfile(cfg.projectRoot,paths(1)),'final_state','-v7');
reread=load(fullfile(cfg.projectRoot,paths(1)),'final_state');assert(isequal(final_state,reread.final_state));
bindings.state_file=paths(1);bindings.state_sha256=q2_validation_io('hash',fullfile(cfg.projectRoot,paths(1)));
bindings.units=struct('time','s','radius','m','temperature','K','moisture','kg/kg');
bindings.time_s=10800;bindings.internal_nodes=numel(solution.radius_m);
bindings.q3_production_restart_authorized=false;
bindings.authorization_rule='Requires effective q2 validation PASS with this exact state binding; stage30 does not grant it.';
q2_validation_io('json',fullfile(cfg.projectRoot,paths(2)),bindings);
end

function paths=plot_q2_results(cfg,fields,mechanism,validation)
p=q2_compute("parameters");[ambient,~]=q2_compute("load_inputs",cfg);times=(0:10800)'/3600;
fieldValues={ [28*ones(1,21);fields.temperature_C], [p.C0*ones(1,21);fields.moisture_kgkg] };
[Ta,Ce]=q2_compute("interpolate_ambient",(0:10800)',ambient);
outside={Ta,Ce};unitLabels=["温度（℃）","干基含水率（kg/kg）"];
stems=["result_q2_01_temperature_field","result_q2_02_moisture_field","result_q2_03_coupling_mechanism"];
paths="figures/q2/result/"+stems+".svg";folder=fullfile(cfg.figuresDir,'q2','result');
if ~isfolder(folder),mkdir(folder);end
qa=cfg.figureQADir;if ~isfolder(qa),mkdir(qa);end
blue=[0,0.4471,0.6980];orange=[0.8353,0.3686,0];gray=[0.3,0.3,0.3];
audits=cell(3,1);
for field=1:2
    fig=figure('Visible','off','Color','w');cleanup=onCleanup(@()close(fig));
    style=q2Style(fig);fig.Units='centimeters';fig.Position(3:4)=[15.5,9.4];
    data=fieldValues{field};lo=min(data,[],'all');hi=max(data,[],'all');
    ax=axes(fig,'Position',[.105,.34,.48,.58]);
    renderRadius=linspace(0,2,800);displayData=interp1(fields.radius_cm,data',renderRadius,'linear');
    imagesc(ax,times,renderRadius,displayData);axis(ax,'xy');clim(ax,[lo,hi]);colormap(ax,blueMap());
    set(ax,'XLim',[0,3],'YLim',[0,2]);decorate(ax,[0,1,2,3],[0,1,2],"时间（小时）","半径（厘米）",style.font);
    label(ax,0,1.06,'(a)',style.font,'left');
    bar=axes(fig,'Position',[.105,.16,.48,.025]);
    imagesc(bar,[lo,hi],[0,1],linspace(lo,hi,1200));axis(bar,'xy');colormap(bar,blueMap());
    set(bar,'XLim',[lo,hi],'YLim',[0,1],'XTick',[],'YTick',[]);
    for v=[lo,(lo+hi)/2,hi]
        text(bar,v,-1.2,sprintf('%.3g',v),'HorizontalAlignment','center','FontName',style.font, ...
            'FontSize',7.5,'Interpreter','none','Clipping','off');
    end
    label(bar,.5,-3.5,unitLabels(field),style.font,'center');
    bar.Box='on';bar.XAxis.Visible='off';bar.YAxis.Visible='off';
    b=axes(fig,'Position',[.715,.34,.245,.58]);hold(b,'on');
    plot(b,times,data(:,1),'-','Color',blue);plot(b,times,data(:,end),'--','Color',orange);
    plot(b,times,outside{field},':','Color',gray);xlim(b,[0,3]);
    if field==1,ylim(b,[27,52]);yt=[28,40,50];else,ylim(b,[0,2.65]);yt=[0,1,2];end
    decorate(b,[0,1.5,3],yt,"时间（小时）",unitLabels(field),style.font);label(b,0,1.06,'(b)',style.font,'left');
    legendAxes=axes(fig,'Position',[.68,.08,.30,.07],'Visible','off');hold(legendAxes,'on');
    manualLegend(legendAxes,[blue;orange;gray],{'-','--',':'},["中心","表面","环境"],style.font);
    drawnow;assert(isempty(audit_publication_figure(fig)),'CUMCM:q2:FigureLayout','Figure layout check failed.');
    export_publication_figure(fig,fullfile(folder,stems(field)),Formats="svg",DPI=300,QADir=qa, ...
        StrictDesign=true,Width=15.5,Height=9.4,Units="centimeters");
    audits{field}=q2_editable_svg(fig,fullfile(cfg.projectRoot,paths(field)),15.5,9.4);
    clear cleanup
end
fig=figure('Visible','off','Color','w');cleanup=onCleanup(@()close(fig));
style=q2Style(fig);fig.Units='centimeters';fig.Position(3:4)=[15.5,11.5];
initial=q2_compute("material_properties",p.C0,p.T0_K,p);
a=axes(fig,'Position',[.12,.55,.83,.34]);hold(a,'on');
plot(a,times,[initial.D;mechanism.D(:,1)]*1e9,'-','Color',blue);
plot(a,times,[initial.D;mechanism.D(:,2)]*1e9,'--','Color',orange);xlim(a,[0,3]);
decorate(a,[0,1,2,3],a.YTick,"时间（小时）","水分扩散系数（10^-9 m²/s）",style.font);label(a,0,1.08,'(a)',style.font,'left');
lg=axes(fig,'Position',[.62,.91,.3,.05],'Visible','off');hold(lg,'on');
manualLegend(lg,[blue;orange],{'-','--'},["中心","表面"],style.font);
b=axes(fig,'Position',[.12,.13,.34,.26]);hold(b,'on');
plot(b,times,[initial.alpha;mechanism.alpha(:,1)]*1e7,'-','Color',blue);
plot(b,times,[initial.alpha;mechanism.alpha(:,2)]*1e7,'--','Color',orange);xlim(b,[0,3]);
decorate(b,[0,1.5,3],b.YTick,"时间（小时）","热扩散率（10^-7 m²/s）",style.font);label(b,0,1.12,'(b)',style.font,'left');
c=axes(fig,'Position',[.63,.13,.32,.26]);hold(c,'on');
plot(c,times,[0;mechanism.logD_C(:,2)],'--','Color',orange);
plot(c,times,[0;mechanism.logD_T(:,2)],'-','Color',blue);
plot(c,times,[0;mechanism.logD_C(:,2)+mechanism.logD_T(:,2)],':','Color',gray);xlim(c,[0,3]);
decorate(c,[0,1.5,3],c.YTick,"时间（小时）","对数相对变化（无量纲）",style.font);label(c,0,1.12,'(c)',style.font,'left');
lg2=axes(fig,'Position',[.59,.445,.38,.045],'Visible','off');hold(lg2,'on');
manualLegend(lg2,[blue;orange;gray],{'-','--',':'},["升温贡献","失水贡献","合计"],style.font);
drawnow;assert(isempty(audit_publication_figure(fig)),'CUMCM:q2:FigureLayout','Figure layout check failed.');
export_publication_figure(fig,fullfile(folder,stems(3)),Formats="svg",DPI=300,QADir=qa, ...
    StrictDesign=true,Width=15.5,Height=11.5,Units="centimeters");
audits{3}=q2_editable_svg(fig,fullfile(cfg.projectRoot,paths(3)),15.5,11.5);
q2_validation_io('json',fullfile(cfg.projectRoot,'artifacts','q2','figure_export_audit.json'), ...
    struct('run_id',validation.run_id,'figures',{audits},'status','CANDIDATE_AWAITING_VISUAL_REVIEW', ...
    'data_rows_per_field',10800,'radial_columns',21,'time_values_discarded',0, ...
    'radial_display_resampling','21 samples linearly interpolated to800 pixels; output fields unchanged', ...
    'formal_document_insertion','NOT_RUN_STAGE30_CANDIDATE_DIMENSIONS'));
end
function decorate(ax,xt,yt,xname,yname,font)
ax.XTick=xt;ax.YTick=yt;ax.XTickLabel={};ax.YTickLabel={};ax.Box='off';
ax.XAxis.Exponent=0;ax.YAxis.Exponent=0;ax.TickDir='out';ax.FontName=font;
for x=xt
    text(ax,x,ax.YLim(1)-.075*diff(ax.YLim),sprintf('%.3g',x),'HorizontalAlignment','center', ...
        'VerticalAlignment','top','FontName',font,'FontSize',7.5,'Interpreter','none','Clipping','off');
end
for y=yt
    text(ax,ax.XLim(1)-.035*diff(ax.XLim),y,sprintf('%.3g',y),'HorizontalAlignment','right', ...
        'VerticalAlignment','middle','FontName',font,'FontSize',7.5,'Interpreter','none','Clipping','off');
end
if strlength(xname)>0,label(ax,.5,-.22,xname,font,'center');end
if strlength(yname)>0
    text(ax,-.058/ax.Position(3),.5,yname,'Units','normalized','Rotation',90,'HorizontalAlignment','center', ...
        'FontName',font,'FontSize',8,'Interpreter','none','Clipping','off');
end
end
function label(ax,x,y,value,font,alignment)
text(ax,x,y,value,'Units','normalized','FontName',font,'FontSize',8,'Interpreter','none', ...
    'HorizontalAlignment',alignment,'Clipping','off');
end
function manualLegend(ax,colors,styles,names,font)
xlim(ax,[0,1]);ylim(ax,[0,1]);ax.Visible='off';
for k=1:numel(names)
    x=(k-1)/numel(names);plot(ax,x+[0,.07],[.5,.5],styles{k},'Color',colors(k,:),'LineWidth',1.1);
    text(ax,x+.085,.5,names(k),'FontName',font,'FontSize',7.5,'Interpreter','none','VerticalAlignment','middle');
end
end
function style=q2Style(fig)
style=apply_publication_style(fig,'en','report');
assert(any(string(listfonts)=="Microsoft YaHei UI"),'CUMCM:q2:FigureFont','Required CJK font unavailable.');
style.font="Microsoft YaHei UI";set(fig,'DefaultAxesFontName',style.font,'DefaultTextFontName',style.font);
end
function colors=blueMap()
colors=interp1([0,1],[.95,.97,.99;.05,.20,.38],linspace(0,1,256));
end

function audit=q2_editable_svg(fig,destination,widthCm,heightCm)
fig.Units='centimeters';fig.Position(3:4)=[widthCm,heightCm];
fig.PaperUnits='centimeters';fig.PaperPosition=[0,0,widthCm,heightCm];
drawnow;texts=findall(fig,'Type','text','Visible','on');
keep=arrayfun(@(t)strlength(join(string(t.String)))>0,texts);texts=texts(keep);
fig.Units='pixels';sizePx=fig.Position(3:4);fig.Units='centimeters';
records=struct([]);
for k=1:numel(texts)
    t=texts(k);assert(string(t.Interpreter)=="none",'CUMCM:q2:SvgText','Editable adapter requires plain text.');
    old=t.Units;t.Units='pixels';pos=t.Position;t.Units=old;
    ax=ancestor(t,'axes');assert(~isempty(ax),'CUMCM:q2:SvgText','Unsupported text parent.');
    box=getpixelposition(ax,true);record=struct('x',(box(1)+pos(1)-1)/sizePx(1), ...
        'y',1-(box(2)+pos(2)-1)/sizePx(2),'text',string(t.String), ...
        'font',string(t.FontName),'size',t.FontSize,'weight',string(t.FontWeight), ...
        'color',t.Color,'rotation',t.Rotation,'align',string(t.HorizontalAlignment), ...
        'valign',string(t.VerticalAlignment));
    if isempty(records),records=record;else,records(k)=record;end
    t.Visible='off';
end
restore=onCleanup(@()set(texts,'Visible','on'));
native=string(tempname(fileparts(destination)))+".svg";print(fig,char(native),'-dsvg','-painters','-r300');
xml=fileread(native);
imageObjects=flip(findall(fig,'Type','image'));
tokens=regexp(xml,'(<image\s[^>]*?(?:xlink:)?href=")([^"]*)(")','tokens');
assert(numel(tokens)==numel(imageObjects),'CUMCM:q2:SvgImage','Image-layer count differs.');
pixelSizes=zeros(numel(imageObjects),2);
for k=1:numel(imageObjects)
    im=imageObjects(k);ax=ancestor(im,'axes');data=im.CData;
    if ndims(data)==2
        colors=colormap(ax);index=1+floor((double(data)-ax.CLim(1))/diff(ax.CLim)*size(colors,1));
        index=min(size(colors,1),max(1,index));rgb=uint8(round(255*reshape(colors(index(:),:),[size(data),3])));
    else
        rgb=data;
    end
    if string(ax.YDir)=="normal",rgb=flipud(rgb);end
    if size(rgb,1)<64,rgb=repelem(rgb,ceil(64/size(rgb,1)),1,1);end
    png=string(tempname(fileparts(destination)))+".png";imwrite(rgb,png);
    fid=fopen(png,'rb');bytes=fread(fid,Inf,'*uint8');fclose(fid);
    encoded=char(matlab.net.base64encode(bytes));
    old=strjoin(tokens{k},'');replacement=[tokens{k}{1},'data:image/png;base64,',encoded,tokens{k}{3}];
    xml=replace(xml,old,replacement);pixelSizes(k,:)=[size(rgb,2),size(rgb,1)];
    internal=fullfile(fileparts(fileparts(fileparts(fileparts(destination)))),'artifacts','q2','svg_native');
    if ~isfolder(internal),mkdir(internal);end
    [~,stemImage]=fileparts(destination);movefile(png,fullfile(internal,string(stemImage)+"_layer"+k+".png"),'f');
end
w=regexp(xml,'<svg[^>]*\swidth="([\d.]+)"','tokens','once');
h=regexp(xml,'<svg[^>]*\sheight="([\d.]+)"','tokens','once');
assert(~isempty(w)&&~isempty(h),'CUMCM:q2:SvgText','SVG coordinate dimensions missing.');
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
qaFolder=fullfile(fileparts(fileparts(fileparts(fileparts(destination)))),'artifacts','q2','svg_native');
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
    'width_cm',widthCm,'height_cm',heightCm,'native_matlab_data_geometry_preserved',true, ...
    'embedded_image_pixels_width_height',pixelSizes,'preview_from_final_svg',true);
end
function s=escape(s)
s=replace(string(s),["&","<",">",'"'],["&amp;","&lt;","&gt;","&quot;"]);
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
            " 个坐标轴的标签或刻度可能超出画布";
    end

    titleParts = string(ax.Title.String);
    titleText = strtrim(join(titleParts(:), " "));
    if strlength(titleText) > 30
        issues(end+1) = "第 " + k + ...
            " 个坐标轴标题超过 30 字符；只允许必要的短面板标题，完整说明写入论文题注";
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
                " 个坐标轴对稠密曲线逐点绘制标记；应取消标记或设置 MarkerIndices";
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
                " 个坐标轴的柱状图未从零开始；若必须截断，应改用点图/区间图并显式说明";
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
                " 个坐标轴是已标数值的 2×2 小矩阵，不应再使用冗余 colorbar";
            break;
        end
    end
end

legendList = findall(fig, "Type", "legend");
for legendObject = reshape(legendList, 1, [])
    labels = string(legendObject.String);
    if numel(labels) > 5
        issues(end+1) = "图例超过 5 项；应直接标注、改为共享图例或拆分证据";
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
    originalPaperUnits, originalPaperPosition));
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
    formalPaths(end+1,1) = pathValue;
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
