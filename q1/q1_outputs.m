function varargout = q1_outputs(action,varargin)
action=string(action);
assert(isscalar(action),"CUMCM:q1:InvalidOutputAction","输出动作必须是标量字符串。");
switch action
    case "write_csv"
        varargout{1}=write_q1_csv(varargin{:});
    case "write_result1"
        varargout{1}=write_q1_result1(varargin{:});
    case "plot_results"
        varargout{1}=plot_q1_results(varargin{:});
    otherwise
        error("CUMCM:q1:InvalidOutputAction","未知输出动作：%s",action);
end
end

function paths=write_q1_csv(cfg,fields,validation)
folder=fullfile(cfg.resultsDir,'q1'); if ~isfolder(folder), mkdir(folder); end
names=["q1_temperature_unrounded.csv","q1_moisture_unrounded.csv", ...
    "q1_table1_temperature.csv","q1_table2_moisture.csv"];
data={fields.temperature_C,fields.moisture_kgkg, ...
    fields.temperature_C(fields.paper_time_s,fields.paper_columns), ...
    fields.moisture_kgkg(fields.paper_time_s,fields.paper_columns)};
paths=strings(1,4);
for k=1:4
    if k<=2, time=fields.time_s; radius=fields.radius_cm;
    else, time=fields.paper_time_s; radius=fields.radius_cm(fields.paper_columns); end
    if k<=2, matrix=data{k}; format=',%.17g';
    else, matrix=round(data{k},4); format=',%.4f'; end
    path=fullfile(folder,names(k)); fid=fopen(path,'w','n','UTF-8'); assert(fid>=0);
    clean=onCleanup(@()fclose(fid));
    fprintf(fid,'%s','时间\到药材中心的距离'); fprintf(fid,',%.1f',radius); fprintf(fid,'\n');
    for row=1:numel(time)
        fprintf(fid,'%.17g',time(row)); fprintf(fid,format,matrix(row,:)); fprintf(fid,'\n');
    end
    clear clean
    actual=readmatrix(path,'NumHeaderLines',1);
    assert(isequal(actual,[time(:),matrix]),'CUMCM:q1:OutputMappingFailure','CSV round-trip mismatch.');
    paths(k)="results/q1/"+names(k);
end
fid=fopen(fullfile(folder,'q1_validation_summary.json'),'w','n','UTF-8'); assert(fid>=0);
clean=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(validation,'PrettyPrint',true));
paths(end+1)="results/q1/q1_validation_summary.json";
end

function path=write_q1_result1(cfg,fields)
source=fullfile(cfg.rawDataDir,'附件3','result1.xlsx'); originalHash=hashFile(source);
folder=fullfile(cfg.processedDataDir,'q1'); if ~isfolder(folder), mkdir(folder); end
work=fullfile(cfg.projectRoot,'artifacts','q1','output_work'); if ~isfolder(work), mkdir(work); end
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
    rows=cell(1801,1); header=cell(1,22);
    header{1}='<c r="A1" s="2" t="s"><v>1</v></c>';
    headStyle=3-sheet; dataStyle=3+sheet;
    for col=2:22
        header{col}=sprintf('<c r="%s1" s="%d"><v>%.1f</v></c>',char(64+col),headStyle,fields.radius_cm(col-1));
    end
    rows{1}=['<row r="1" spans="1:22" ht="14.1" customHeight="1">',strjoin(header,''),'</row>'];
    for time=1:1800
        cells=cell(1,22); row=time+1;
        cells{1}=sprintf('<c r="A%d" s="1"><v>%d</v></c>',row,time);
        for col=2:22
            cells{col}=sprintf('<c r="%s%d" s="%d"><v>%.4f</v></c>',char(64+col),row,dataStyle,values{sheet}(time,col-1));
        end
        rows{row}=sprintf('<row r="%d" spans="1:22" ht="14.1" customHeight="1">%s</row>',row,strjoin(cells,''));
    end
    xml=regexprep(xml,'<sheetData>.*?</sheetData>',['<sheetData>',strjoin(rows,''),'</sheetData>']);
    xml=replace(xml,'ref="A1:F5"','ref="A1:V1801"');
    xml=replace(xml,'min="2" max="6"','min="2" max="22"');
    xml=replace(xml,'min="7" max="16384"','min="23" max="16384"');
    writeText(pathXml,xml);
end
shared=fullfile(parts,'xl','sharedStrings.xml'); ss=fileread(shared);
ss=regexprep(ss,' count="\d+"',' count="2"','once'); writeText(shared,ss);
first=string(tempname(work))+".xlsx"; second=string(tempname(work))+".xlsx";
q1_pack_ooxml(parts,first); canonicalize_q1_xlsx(first,cfg.projectRoot);
q1_pack_ooxml(parts,second); canonicalize_q1_xlsx(second,cfg.projectRoot);
assert(strcmp(hashFile(first),hashFile(second)),'CUMCM:q1:NonDeterministicArtifact','Workbook writes differ.');
sheets=sheetnames(first); assert(isequal(sheets(:),["温度";"水分浓度"]));
for k=1:2
    top=readcell(first,'Sheet',sheets(k),'Range','A1','UseExcel',false);
    assert(string(top{1})=="时间\到药材中心的距离");
    matrix=readmatrix(first,'Sheet',sheets(k),'Range','A1:V1801','UseExcel',false);
    assert(isequal(matrix(1,2:end),fields.radius_cm) && isequal(matrix(2:end,1),fields.time_s));
    assert(isequal(matrix(2:end,2:end),values{k}),'CUMCM:q1:OutputMappingFailure','XLSX values mismatch.');
end
assert(strcmp(hashFile(source),originalHash),'CUMCM:q1:TemplateMutation','Raw template changed.');
output=fullfile(folder,'result1.xlsx');
if isfile(output), copyfile(output,string(tempname(work))+"_previous.xlsx"); end
movefile(first,output,'f'); path="data/processed/q1/result1.xlsx";
end
function writeText(path,text)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0); c=onCleanup(@()fclose(fid)); fprintf(fid,'%s',text);
end
function value=hashFile(path)
bytes=java.nio.file.Files.readAllBytes(java.io.File(char(path)).toPath());
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end

function canonicalize_q1_xlsx(path,projectRoot)
root=string(java.io.File(char(projectRoot)).getCanonicalPath());
actual=string(java.io.File(char(path)).getCanonicalPath());
assert(startsWith(lower(actual),lower(root+filesep)) && ~contains(lower(actual),lower(filesep+"data"+filesep+"raw"+filesep)), ...
    'CUMCM:q1:TemplateMutation','Only project-generated workbooks may be canonicalized.');
work=fullfile(root,'artifacts','q1','output_work'); if ~isfolder(work), mkdir(work); end
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
q1_pack_ooxml(parts,path);
end

function q1_pack_ooxml(parts,destination)
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

function paths=plot_q1_results(cfg,fields,validation)
p=q1_compute("parameters"); [ambient,~]=q1_compute("load_inputs",cfg);
time=(0:1800)'/60; radius=fields.radius_cm;
data={fields.temperature_C,fields.moisture_kgkg}; initials=[p.theta0,p.C0];
boundary={interp1(ambient.time_s,ambient.temperature_C,(0:1800)'), ...
    interp1(ambient.time_s,ambient.moisture_kgkg,(0:1800)')};
names=["result_q1_temperature_field","result_q1_moisture_field"];
labels=["温度（℃）","干基含水率（kg/kg）"];
paths=strings(1,2); audits=cell(1,2);
for k=1:2
    assert(isequal(size(data{k}),[1800,21])&&all(isfinite(data{k}),'all'));
    fig=figure('Visible','off','Color','white'); closeFig=onCleanup(@()close(fig));
    style=apply_publication_style(fig,"zh","report");
    fig.Units='centimeters'; fig.Position(3:4)=[15.5,9.4];
    ax=axes(fig,'Position',[0.12,0.51,0.69,0.42]);
    matrix=[initials(k)*ones(1,21);data{k}];
    contourf(ax,time,radius,matrix',24,'LineColor','none');
    anchors=[68,1,84;71,44,122;59,81,139;44,113,142;33,144,141; ...
        39,173,129;92,200,99;170,220,50;253,231,37]/255;
    colormap(ax,interp1(linspace(0,1,9),anchors,linspace(0,1,256)));
    set(ax,'YDir','normal','Box','off');
    xlim(ax,[0,30]); ylim(ax,[0,2]); xticks(ax,0:5:30); yticks(ax,0:0.5:2);
    xlabel(ax,'时间（min）','Interpreter','none'); ylabel(ax,'半径（cm）','Interpreter','none');
    cb=colorbar(ax); cb.Position=[0.855,0.51,0.022,0.42];
    cb.Label.String=labels(k); cb.Label.Interpreter='none'; cb.FontSize=7.5;
    cb.FontName=style.font; cb.Label.FontName=style.font;
    ax.Position=[0.12,0.51,0.69,0.42];
    small=axes(fig,'Position',[0.12,0.13,0.69,0.22]); hold(small,'on');
    plot(small,time,matrix(:,1),'-','Color',style.colors(1,:),'LineWidth',1.1);
    plot(small,time,matrix(:,end),'--','Color',style.colors(4,:),'LineWidth',1.1);
    plot(small,time,boundary{k},':','Color',style.colors(3,:),'LineWidth',1.25);
    xlabel(small,'时间（min）','Interpreter','none'); ylabel(small,labels(k),'Interpreter','none');
    xlim(small,[0,30]); xticks(small,0:5:30); grid(small,'off'); small.Box='off';
    if k==1, ylim(small,[27.5,42]); external='环境温度';
    else, ylim(small,[0,2.75]); external='等效外界'; end
    lg=legend(small,{'中心','表面',external},'Orientation','horizontal','Box','off');
    lg.FontName=style.font; lg.FontSize=7.5; lg.Units='normalized'; lg.Position=[0.20,0.375,0.57,0.06];
    stem=fullfile(cfg.figuresDir,'q1','result',names(k));
    audits{k}=export_publication_figure(fig,stem,Formats="svg",DPI=300, ...
        QADir=cfg.figureQADir,Width=15.5,Height=9.4,Units="centimeters",StrictDesign=true);
    paths(k)="figures/q1/result/"+names(k)+".svg";
    clear closeFig
end
profile=struct('kind','DRAFT_FIGURE_PROFILE','status',validation.status, ...
    'rows_per_field',1800,'radius_columns',21,'excluded_points',0,'all_finite',true, ...
    'temperature_range_C',[min(data{1},[],'all'),max(data{1},[],'all')], ...
    'moisture_range_kgkg',[min(data{2},[],'all'),max(data{2},[],'all')], ...
    'width_cm',15.5,'height_cm',9.4,'n_logical_figures',2, ...
    'visual_interpolation','Filled contour visualization of all output nodes; numerical tables are not interpolated.', ...
    'initial_row','Known specified initial state is included at t=0 in figures only.');
folder=fullfile(cfg.projectRoot,'artifacts','q1');
fid=fopen(fullfile(folder,'q1_figure_profile.json'),'w','n','UTF-8'); assert(fid>=0);
c=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(profile,'PrettyPrint',true));
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
