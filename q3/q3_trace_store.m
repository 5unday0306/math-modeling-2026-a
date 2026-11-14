function trace=q3_trace_store(action,trace,payload)
switch char(action)
    case 'init'
        assert(~isfile(trace.path),'CUMCM:q3:TraceExists','Refusing to overwrite a trace.');
        folder=fileparts(trace.path);if ~isfolder(folder),mkdir(folder);end
        trace.width=numel(payload);trace.written=0;trace.used=0;
        trace.buffer=zeros(trace.width,128);
        h5create(trace.path,'/accepted',[trace.width,Inf],'ChunkSize',[trace.width,128],'Datatype','double');
        h5writeatt(trace.path,'/accepted','complete',0);
        trace=q3_trace_store('append',trace,payload);
    case 'append'
        assert(numel(payload)==trace.width&&all(isfinite(payload)));
        assert(8*trace.width*(trace.written+trace.used+1)+1024^2<=trace.maxBytes, ...
            'CUMCM:q3:ResourceLimit','Trace exceeds budget.');
        trace.used=trace.used+1;trace.buffer(:,trace.used)=payload(:);
        if trace.used==size(trace.buffer,2),trace=flush(trace);end
    case {'checkpoint','finalize'}
        trace=flush(trace);
        h5writeatt(trace.path,'/accepted','columns_written',trace.written);
        h5writeatt(trace.path,'/accepted','complete',double(strcmp(action,'finalize')));
    otherwise
        error('CUMCM:q3:TraceAction','Unknown trace action.');
end
end
function trace=flush(trace)
if trace.used>0
    h5write(trace.path,'/accepted',trace.buffer(:,1:trace.used), ...
        [1,trace.written+1],[trace.width,trace.used]);
    trace.written=trace.written+trace.used;trace.used=0;
end
end
