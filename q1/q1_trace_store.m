classdef q1_trace_store < handle
    properties
        file
        buffer
        used = 0
        written = 0
        width
        maxRows
        chunkRows
        otherBytes
        traceLimit
    end
    methods
        function obj=q1_trace_store(file,width,maxRows,p)
            obj.file=string(file); obj.width=width; obj.maxRows=maxRows;
            obj.chunkRows=min(256,maxRows);
            obj.traceLimit=p.maxTraceBytes;
            assert(2*8*width*obj.chunkRows < p.maxWorkingBytes, ...
                "CUMCM:q1:ResourceLimit","Trace buffer exceeds working budget.");
            folder=fileparts(obj.file);
            if ~isfolder(folder), mkdir(folder); end
            assert(~isfile(obj.file),"CUMCM:q1:ResourceLimit","Trace path collision.");
            entries=dir(fullfile(folder,"*.h5"));
            obj.otherBytes=sum([entries.bytes]);
            assert(obj.otherBytes+8*width*maxRows+1024^2<=p.maxTraceBytes, ...
                "CUMCM:q1:ResourceLimit","Required full-precision traces exceed 8 GiB.");
            obj.buffer=zeros(width,obj.chunkRows);
            h5create(obj.file,"/trajectory",[width,Inf],"ChunkSize",[width,obj.chunkRows], ...
                "Datatype","double");
        end
        function append(obj,row)
            assert(numel(row)==obj.width && all(isfinite(row)), ...
                "CUMCM:q1:RangeOrFiniteFailure","Invalid trace row.");
            assert(obj.written+obj.used<obj.maxRows, ...
                "CUMCM:q1:ResourceLimit","Trace step count exceeded planned storage.");
            obj.used=obj.used+1; obj.buffer(:,obj.used)=row(:);
            if obj.used==obj.chunkRows, obj.flush(); end
        end
        function flush(obj)
            if obj.used==0, return; end
            h5write(obj.file,"/trajectory",obj.buffer(:,1:obj.used), ...
                [1,obj.written+1],[obj.width,obj.used]);
            obj.written=obj.written+obj.used; obj.used=0;
            info=dir(obj.file);
            assert(obj.otherBytes+info.bytes<=obj.traceLimit, ...
                "CUMCM:q1:ResourceLimit","Actual trace bytes exceed budget.");
        end
        function finish(obj)
            obj.flush();
        end
    end
end
