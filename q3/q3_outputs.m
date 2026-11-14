function varargout=q3_outputs(action,varargin)
action=string(action);
assert(isscalar(action),'CUMCM:q3:InvalidOutputAction','Invalid output action.');
switch action
    case "finish"
        varargout{1}=finish_q3_weighted_outputs(varargin{:});
    case "write_outputs"
        varargout{1}=write_weighted_outputs(varargin{:});
    case "write_result3"
        varargout{1}=write_q3_result3(varargin{:});
    case "plot_results"
        varargout{1}=plot_q3_results(varargin{:});
    case "plot_selection"
        varargout{1}=plot_weighted_selection(varargin{:});
    otherwise
        error('CUMCM:q3:InvalidOutputAction','Unknown output action: %s',action);
end
end

function result=finish_q3_weighted_outputs(cfg)
if nargin<1
    root=string(fileparts(mfilename('fullpath')));
    cfg=q3_compute("config",root);
end
root=cfg.projectRoot;
loaded=load(fullfile(root,'附件','未舍入结果.mat'));
candidate=loaded.candidate;selection=loaded.selection;
paths=write_weighted_outputs(cfg,candidate,selection,loaded.history);
view=struct('maximum_history',loaded.history,'temporal_difference_s',[], ...
    'figure_note','温度与水分分别加权选窗（0.8/0.2）；热湿耦合数值模拟。');
imagePaths=plot_q3_results(cfg,candidate.fields,candidate.event,view,[1,2]);
selectionPaths=plot_weighted_selection(cfg,loaded.ambient,selection);
pngs=[imagePaths(endsWith(imagePaths,'.png')),selectionPaths(endsWith(selectionPaths,'.png'))];
for k=1:numel(pngs)
    [~,name,extension]=fileparts(pngs(k));
    copyfile(fullfile(root,pngs(k)),fullfile(root,'图片',name+extension));
end
result=struct('method','weighted_variance_and_sample_count', ...
    'variance_weight',selection.temperature.alpha,'sample_weight',1-selection.temperature.alpha, ...
    'temperature_start_min',selection.temperature.start_min,'temperature_samples',selection.temperature.sample_count, ...
    'moisture_start_min',selection.moisture.start_min,'moisture_samples',selection.moisture.sample_count, ...
    'temperature_boundary_C',candidate.scenario.temperature_C,'moisture_boundary_kgkg',candidate.scenario.moisture_kgkg, ...
    'stop_time_h',candidate.event.stop_time_h,'stop_time_s',candidate.event.stop_time_s, ...
    'maximum_moisture_at_stop',candidate.event.maximum_C,'output_rows',numel(candidate.fields.time_s), ...
    'compute_seconds',loaded.computeSeconds,'files',[paths,imagePaths,selectionPaths]);
fid=fopen(fullfile(root,'附件','计算摘要.json'),'w','n','UTF-8');assert(fid>=0);
writer=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,'PrettyPrint',true));
end

function paths=write_weighted_outputs(cfg,candidate,selection,history)
root=cfg.projectRoot;folder=fullfile(root,'附件','q3');if ~isfolder(folder),mkdir(folder);end
fields=candidate.fields;event=candidate.event;
assert(all(isfinite(fields.moisture_kgkg),'all')&&event.maximum_C<.15);
paths=["附件/q3/完整含水率_未舍入.csv","附件/q3/表5_每6小时.csv", ...
    "附件/温度候选窗口.csv","附件/水分浓度候选窗口.csv", ...
    "附件/q3/全空间最大值.csv","附件/边界与结果摘要.csv"];
header="时间\到药材中心的距离";
for r=fields.radius_cm,header=header+","+sprintf('%.1f',r);end
numericCsv(fullfile(root,paths(1)),header,[fields.time_s,fields.moisture_kgkg],'%.17g');
numericCsv(fullfile(root,paths(2)), ...
    '时间（小时）,中心（0厘米）,0.5厘米,1厘米,1.5厘米,表面（2厘米）', ...
    [fields.paper_time_h,fields.paper_moisture_kgkg],'%.4f');
candidateHeader='窗口起点（分钟）,样本数,样本方差,算术平均值,时间平均值,方差惩罚,样本量惩罚,总损失';
numericCsv(fullfile(root,paths(3)),candidateHeader,selection.temperature.candidate_statistics,'%.17g');
numericCsv(fullfile(root,paths(4)),candidateHeader,selection.moisture.candidate_statistics,'%.17g');
numericCsv(fullfile(root,paths(5)),'时间（秒）,全空间最大含水率（kg/kg）',history,'%.17g');
numericCsv(fullfile(root,paths(6)), ...
    '方差权重,样本量权重,温度窗口起点（分钟）,水分窗口起点（分钟）,温度边界（摄氏度）,水分浓度边界（kg/kg）,停止时间（小时）,停止时间（秒）,停止时最大含水率（kg/kg）', ...
    [selection.temperature.alpha,1-selection.temperature.alpha,selection.temperature.start_min,selection.moisture.start_min, ...
    candidate.scenario.temperature_C,candidate.scenario.moisture_kgkg,event.stop_time_h,event.stop_time_s,event.maximum_C],'%.17g');
paths(end+1)=write_q3_result3(cfg,fields);
end
function numericCsv(file,header,data,format)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',header);pattern=[format,repmat([',' format],1,size(data,2)-1),'\n'];fprintf(fid,pattern,data');clear cleanup
if strcmp(format,'%.17g'),assert(isequal(readmatrix(file,'NumHeaderLines',1),data));end
end

function path=write_q3_result3(cfg,fields)
source=fullfile(cfg.rawDataDir,'附件3','result3.xlsx');original=hashFile(source);
folder=fullfile(cfg.processedDataDir,'q3');if ~isfolder(folder),mkdir(folder);end
work=fullfile(cfg.projectRoot,'artifacts','q3','output_work');if ~isfolder(work),mkdir(work);end
parts=string(tempname(work));mkdir(parts);copied=string(tempname(work))+".xlsx";
copyfile(source,copied);unzip(copied,parts);
count=numel(fields.time_s);assert(count+1<=1048576);
stylesPath=fullfile(parts,'xl','styles.xml');styles=fileread(stylesPath);
token=regexp(styles,'<cellXfs count="4">(.*?)</cellXfs>','tokens','once');assert(~isempty(token));
xfs=regexp(token{1},'<xf\s[^>]*>.*?</xf>','match');assert(numel(xfs)==4);
numericStyle=replace(xfs{2},'numFmtId="0"','numFmtId="164"');
numericStyle=replace(numericStyle,'<xf ','<xf applyNumberFormat="1" ');
styles=regexprep(styles,'(<styleSheet[^>]*>)', ...
    '$1<numFmts count="1"><numFmt numFmtId="164" formatCode="0.0000"/></numFmts>','once');
styles=replace(styles,['<cellXfs count="4">',token{1},'</cellXfs>'], ...
    ['<cellXfs count="5">',token{1},numericStyle,'</cellXfs>']);writeText(stylesPath,styles);
xmlPath=fullfile(parts,'xl','worksheets','sheet1.xml');xml=fileread(xmlPath);
assert(~contains(xml,'<mergeCells')&&~contains(xml,'<f>'));
rows=cell(count+1,1);header=cell(1,22);
header{1}='<c r="A1" s="1" t="s"><v>2</v></c>';
for col=2:22,header{col}=sprintf('<c r="%s1" s="1"><v>%.1f</v></c>',char(64+col),fields.radius_cm(col-1));end
rows{1}=['<row r="1" spans="1:22" ht="14.1" customHeight="1">',strjoin(header,''),'</row>'];
values=round(fields.moisture_kgkg,4);
for k=1:count
    cells=cell(1,22);row=k+1;
    cells{1}=sprintf('<c r="A%d" s="1"><v>%.17g</v></c>',row,fields.time_s(k));
    for col=2:22,cells{col}=sprintf('<c r="%s%d" s="4"><v>%.4f</v></c>',char(64+col),row,values(k,col-1));end
    rows{row}=sprintf('<row r="%d" spans="1:22" ht="14.1" customHeight="1">%s</row>',row,strjoin(cells,''));
end
xml=regexprep(xml,'<sheetData>.*?</sheetData>',['<sheetData>',strjoin(rows,''),'</sheetData>']);
xml=replace(xml,'ref="A1:F5"',sprintf('ref="A1:V%d"',count+1));
xml=replace(xml,'min="2" max="6"','min="2" max="22"');
xml=replace(xml,'min="7" max="16384"','min="23" max="16384"');writeText(xmlPath,xml);
shared=fullfile(parts,'xl','sharedStrings.xml');ss=fileread(shared);
ss=regexprep(ss,' count="\d+"',' count="1"','once');writeText(shared,ss);
candidate=string(tempname(work))+".xlsx";q2_pack_ooxml(parts,candidate);
assert(isequal(sheetnames(candidate),"Sheet1"));
top=readcell(candidate,'Sheet','Sheet1','Range','A1','UseExcel',false);
assert(string(top{1})=="时间\到药材中心的距离");
matrix=readmatrix(candidate,'Sheet','Sheet1','Range',sprintf('A1:V%d',count+1),'UseExcel',false);
assert(all(abs(matrix(1,2:end)-fields.radius_cm)<1e-12)&&isequal(matrix(2:end,1),fields.time_s));
assert(isequal(matrix(2:end,2:end),values),'CUMCM:q3:WorkbookMapping','Workbook values differ.');
assert(strcmp(hashFile(source),original),'CUMCM:q3:TemplateMutation','Raw result3 template changed.');
destination=fullfile(folder,'result3.xlsx');
if isfile(destination),copyfile(destination,string(tempname(work))+"_previous.xlsx");end
movefile(candidate,destination,'f');path="附件/q3/result3.xlsx";
end
function writeText(path,value)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',value); %#ok<NASGU>
end
function value=hashFile(path)
bytes=java.nio.file.Files.readAllBytes(java.io.File(char(path)).toPath());
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
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

function paths=plot_q3_results(cfg,fields,event,validation,whichFigures)
if nargin<5,whichFigures=1:3;end
stems=["result_q3_01_threshold_event","result_q3_02_radial_profiles","result_q3_03_boundary_and_numerics"];
folder=fullfile(cfg.figuresDir,'q3','result');pngFolder=fullfile(cfg.resultsDir,'q3','png');
if ~isfolder(folder),mkdir(folder);end
if ~isfolder(pngFolder),mkdir(pngFolder);end
paths=strings(1,0);blue=[0,.4471,.698];orange=[.8353,.3686,0];gray=[.32,.32,.32];
for index=reshape(whichFigures,1,[])
    fig=figure('Visible','off','Color','w');cleanup=onCleanup(@()close(fig));
    apply_publication_style(fig,'en','report');
    font="Microsoft YaHei UI";assert(any(string(listfonts)==font));
    set(fig,'DefaultAxesFontName',font,'DefaultTextFontName',font, ...
        'DefaultTextInterpreter','none','DefaultLegendInterpreter','none');
    height=10;if index==3,height=11.5;end
    fig.Units='centimeters';fig.Position(3:4)=[15.5,height];
    if index==1
        history=validation.maximum_history;
        a=axes(fig,'Position',[.14,.48,.81,.43]);hold(a,'on');
        plot(a,history(:,1)/3600,history(:,2),'Color',blue,'DisplayName','全空间重建最大值');
        yline(a,.15,'--','Color',gray,'DisplayName','题目阈值 0.15');
        xlim(a,[0,event.stop_time_h]);ylim(a,[0,max(history(:,2))*1.05]);
        xlabel(a,'时间（小时）');ylabel(a,'最大含水率（kg/kg）');
        legend(a,'Location','northeast','Box','off');panel(a,'(a) 全程');
        b=axes(fig,'Position',[.14,.15,.81,.20]);hold(b,'on');
        near=abs(history(:,1)-event.nominal_root_s)<=120;
        local=history(near,:);assert(size(local,1)>=2);
        plot(b,local(:,1)-event.nominal_root_s,(local(:,2)-.15)*1e6,'Color',blue);
        yline(b,0,'--','阈值','Color',gray,'FontName',font,'FontSize',7);
        xline(b,event.bracket_left_s-event.nominal_root_s,':','Color',gray);
        xline(b,event.bracket_right_s-event.nominal_root_s,':','Color',gray);
        plot(b,event.stop_time_s-event.nominal_root_s,(event.maximum_C-.15)*1e6, ...
            'o','Color',orange,'MarkerFaceColor',orange,'MarkerSize',4);
        xlim(b,[min(local(:,1))-event.nominal_root_s,max(1,event.stop_time_s-event.nominal_root_s)]);
        xlabel(b,'相对名义临界时间（秒）');ylabel(b,'阈值差（百万分之一 kg/kg）');
        panel(b,'(b) 事件邻域：橙点为停止点');
    elseif index==2
        a=axes(fig,'Position',[.14,.18,.81,.70]);hold(a,'on');
        count=numel(fields.paper_time_h);colors=interp1([0,1],[.55,.72,.82;.02,.20,.36],linspace(0,1,max(count,2)));
        styles={'-','--',':','-.'};handles=gobjects(count,1);labels=cell(1,count);
        [found,indices]=ismember(round(fields.paper_time_h*3600*1000),round(fields.time_s*1000));assert(all(found));
        for k=1:count
            values=fields.moisture_kgkg(indices(k),:);
            handles(k)=plot(a,fields.radius_cm,values,'Color',colors(k,:),'LineStyle',styles{mod(k-1,4)+1});
            if k==count,labels{k}=sprintf('结束 %.4f h',fields.paper_time_h(k));set(handles(k),'Color',orange,'LineWidth',1.5);
            else,labels{k}=sprintf('%.0f h',fields.paper_time_h(k));end
        end
        yline(a,.15,':','题目阈值 0.15','Color',gray,'FontName',font,'FontSize',7,'LabelHorizontalAlignment','right');
        xlim(a,[0,2]);ylim(a,[0,max(fields.paper_moisture_kgkg,[],'all')*1.42]);
        xlabel(a,'距药材中心的距离（厘米）');ylabel(a,'干基含水率（kg/kg）');
        set(a,'XTick',[0,.5,1,1.5,2]);
        legendAxes=axes(fig,'Position',[.18,.705,.74,.165],'Visible','off');hold(legendAxes,'on');
        xlim(legendAxes,[0,1]);ylim(legendAxes,[0,1]);rows=ceil(count/2);
        for k=1:count
            col=floor((k-1)/rows);row=mod(k-1,rows);x=.01+.50*col;y=1-(row+.5)/rows;
            plot(legendAxes,x+[0,.065],[y,y],'Color',handles(k).Color,'LineStyle',handles(k).LineStyle,'LineWidth',1.2);
            text(legendAxes,x+.08,y,labels{k},'FontSize',7,'VerticalAlignment','middle','Color',[.15,.15,.15]);
        end
    else
        a=axes(fig,'Position',[.30,.58,.64,.29]);hold(a,'on');
        entries=validation.scenarios;
        done=[entries.completed];values=[entries(done).stop_time_h];positions=find(done);
        plot(a,values,positions,'o','LineStyle','none','Color',blue,'MarkerFaceColor',blue,'MarkerSize',4);
        set(a,'YTick',1:3,'YTickLabel',{entries.label},'YDir','reverse');ylim(a,[.5,3.5]);
        if ~isempty(values)
            span=max(max(values)-min(values),.01);xlim(a,[min(values)-.25*span,max(values)+.45*span]);
            for k=1:numel(values),text(a,values(k)+.025*span,positions(k)-.16,sprintf('%.4f',values(k)),'FontSize',7);end
        end
        xlabel(a,'条件停止时间（小时）');panel(a,'(a) 长期边界情景');
        if ~all(done)
            for missing=find(~done),text(a,mean(a.XLim),missing,'未完成','FontSize',7,'Color',gray);end
        end
        b=axes(fig,'Position',[.30,.19,.64,.23]);hold(b,'on');
        widths=event.bracket_right_s-event.bracket_left_s;
        names={'名义根括区'};
        if ~isempty(validation.temporal_difference_s)
            widths(end+1)=validation.temporal_difference_s;names{end+1}='同初值时间细化';
        end
        barh(b,1:numel(widths),widths,.42,'FaceColor',gray,'EdgeColor','none');
        set(b,'YTick',1:numel(names),'YTickLabel',names,'YDir','reverse');ylim(b,[.5,numel(names)+.5]);
        xlim(b,[0,max(max(widths)*1.4,.03)]);xlabel(b,'区间宽度或临界时间变化（秒）');panel(b,'(b) 数值检查（非统计误差棒）');
        for k=1:numel(widths),text(b,widths(k)+.025*diff(b.XLim),k,sprintf('%.4g',widths(k)),'FontSize',7,'VerticalAlignment','middle');end
        annotation(fig,'textbox',[.12,.065,.82,.04],'String','空间与上游误差未认证；未计算项不作零值。', ...
            'LineStyle','none','FontName',font,'FontSize',7,'HorizontalAlignment','center');
    end
    ax=findall(fig,'Type','axes');for k=1:numel(ax),set(ax(k),'Box','off','FontName',font,'FontSize',7.5,'TickDir','out');grid(ax(k),'off');end
    note='热湿耦合数值模拟。';if isfield(validation,'figure_note'),note=validation.figure_note;end
    annotation(fig,'textbox',[.10,.015,.85,.04],'String',note, ...
        'LineStyle','none','FontName',font,'FontSize',7,'HorizontalAlignment','center');
    drawnow;
    export_publication_figure(fig,fullfile(folder,stems(index)),Formats="svg",DPI=300, ...
        QADir=cfg.figureQADir,StrictDesign=true,Width=15.5,Height=height,Units="centimeters");
    q3_normalize_svg_size(fullfile(folder,stems(index)+".svg"),15.5,height);
    fig.PaperUnits='centimeters';fig.PaperPosition=[0,0,15.5,height];
    png=fullfile(pngFolder,stems(index)+".png");print(fig,char(png),'-dpng','-r300');
    info=imfinfo(png);assert(info.Width>=floor(15.5/2.54*300)&&info.Height>=floor(height/2.54*300));
    assert(isfield(info,'XResolution')&&info.XResolution*.0254>=299.9);
    paths=[paths,"图片源文件/q3/result/"+stems(index)+".svg","附件/q3/png/"+stems(index)+".png"]; %#ok<AGROW>
    clear cleanup
end
end
function panel(ax,label)
text(ax,0,1.12,label,'Units','normalized','FontSize',8,'Clipping','off','VerticalAlignment','bottom');
end

function paths=plot_weighted_selection(cfg,ambient,selection)
fig=figure('Visible','off','Color','w');cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
apply_publication_style(fig,'en','report');font='Microsoft YaHei UI';
set(fig,'DefaultAxesFontName',font,'DefaultTextFontName',font,'DefaultTextInterpreter','none','DefaultLegendInterpreter','none');
fig.Units='centimeters';fig.Position(3:4)=[15.5,13.5];
blue=[0,.4471,.698];orange=[.8353,.3686,0];gray=[.48,.48,.48];
entries={selection.temperature,selection.moisture};values={ambient.temperature_C,ambient.moisture_kgkg};
names={'温度','水分浓度'};units={'温度（℃）','水分浓度（kg/kg）'};
for c=1:2
    entry=entries{c};stats=entry.candidate_statistics;x=.14+.46*(c-1);
    a=axes(fig,'Position',[x,.57,.35,.29]);hold(a,'on');
    h1=plot(a,stats(:,1),stats(:,8),'Color',blue,'LineWidth',1.3,'DisplayName','总损失');
    h2=plot(a,stats(:,1),stats(:,6),'--','Color',orange,'DisplayName','方差惩罚');
    h3=plot(a,stats(:,1),stats(:,7),':','Color',gray,'DisplayName','样本惩罚');
    plot(a,entry.start_min,entry.loss,'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',4,'HandleVisibility','off');
    xline(a,entry.start_min,':','Color',blue,'HandleVisibility','off');
    xlim(a,[entry.first_crossing_min,max(stats(:,1))+1]);ylim(a,[0,max(stats(:,8))*1.07]);
    xlabel(a,'候选窗口起点（分钟）');ylabel(a,'损失（无量纲）');
    text(a,0,1.12,sprintf('(%c) %s：选中%d分钟起点','a'+c-1,names{c},entry.start_min), ...
        'Units','normalized','FontSize',7.5,'VerticalAlignment','bottom');
    if c==1
        lg=legend(a,[h1,h2,h3],'Orientation','horizontal','Box','off','FontSize',7);
        lg.Units='normalized';lg.Position=[.20,.91,.69,.04];
        a.Position=[x,.57,.35,.29];
    end
    b=axes(fig,'Position',[x,.16,.35,.25]);hold(b,'on');
    use=ambient.time_s/60>=entry.first_crossing_min;
    t=ambient.time_s(use)/60;y=values{c}(use);margin=max(max(y)-min(y),eps)*.16;
    limits=[min(y)-margin,max(y)+margin];ylim(b,limits);xlim(b,[t(1),t(end)]);
    patch(b,[entry.start_min,t(end),t(end),entry.start_min], ...
        [limits(1),limits(1),limits(2),limits(2)],[.91,.94,.97],'EdgeColor','none');
    plot(b,t,y,'Color',blue,'LineWidth',.85);
    yline(b,entry.time_mean,'--','Color',orange,'LineWidth',1.1);
    xlabel(b,'观测时间（分钟）');ylabel(b,units{c});
    text(b,0,1.13,sprintf('(%c) 入选%d—240分钟，共%d点','c'+c-1,entry.start_min,entry.sample_count), ...
        'Units','normalized','FontSize',7.3,'VerticalAlignment','bottom');
    if c==1,number=sprintf('时间均值 %.8f ℃',entry.time_mean);
    else,number=sprintf('时间均值 %.11f kg/kg',entry.time_mean);end
    text(b,.5,-.35,number,'Units','normalized','HorizontalAlignment','center','FontSize',7);
    if c==2,b.YAxis.Exponent=0;ytickformat(b,'%.5f');end
end
ax=findall(fig,'Type','axes');for k=1:numel(ax),set(ax(k),'Box','off','FontSize',7,'TickDir','out','FontName',font);end
annotation(fig,'textbox',[.10,.012,.85,.035],'String','权重0.8/0.2；下图蓝线为观测，橙虚线为均值，浅色为入选窗口（非置信区间）。', ...
    'LineStyle','none','FontName',font,'FontSize',6.5,'HorizontalAlignment','center');
stem="result_q3_03_weighted_selection";folder=fullfile(cfg.figuresDir,'q3','result');
export_publication_figure(fig,fullfile(folder,stem),Formats="svg",DPI=300,QADir=cfg.figureQADir, ...
    StrictDesign=true,Width=15.5,Height=13.5,Units="centimeters");
q3_normalize_svg_size(fullfile(folder,stem+".svg"),15.5,13.5);
pngFolder=fullfile(cfg.resultsDir,'q3','png');if ~isfolder(pngFolder),mkdir(pngFolder);end
fig.PaperUnits='centimeters';fig.PaperPosition=[0,0,15.5,13.5];print(fig,fullfile(pngFolder,stem),'-dpng','-r300');
paths=["图片源文件/q3/result/"+stem+".svg","附件/q3/png/"+stem+".png"];
end

function q3_normalize_svg_size(file,width_cm,height_cm)
text=fileread(file);header=regexp(text,'<svg\s[^>]*>','match','once');assert(~isempty(header));
updated=header;
if ~contains(updated,'viewBox=')
    canvas=regexp(text,'<path d="M0 0 L([0-9.]+) 0 L[0-9.]+ ([0-9.]+) L0 ','tokens','once');
    assert(numel(canvas)==2,'CUMCM:q3:SvgCanvas','Cannot identify original SVG canvas.');
    updated=regexprep(updated,'>$',sprintf(' viewBox="0 0 %s %s">',canvas{1},canvas{2}));
end
updated=regexprep(updated,'\swidth="[^"]*"',sprintf(' width="%.15gcm"',width_cm),'once');
updated=regexprep(updated,'\sheight="[^"]*"',sprintf(' height="%.15gcm"',height_cm),'once');
assert(~strcmp(header,updated)||contains(header,sprintf('width="%.15gcm"',width_cm)));
text=replace(text,header,updated);fid=fopen(file,'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',text); %#ok<NASGU>
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
