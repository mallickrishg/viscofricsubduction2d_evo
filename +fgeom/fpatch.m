classdef fpatch < handle
    
    properties
        N;
        % identification number
        id;
        % patch position (upper corner and centre)
        x;
        xc;
        % patch dimension (length and width)
        W;
        % patch orientation (degrees)
        dip;
        % patch slip (amplitude)
        slip;
        % unit vectors in strike, dip, normal and rake directions
        dv;nv;        
        % earth model
        earthModel;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = fpatch()            
            
            if (0==nargin)
                return
            end
            
        end % constructor                                               
        
    end % methods
    
    methods(Static)
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %   convert segment definition to fault patches    %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function [fm] = seg2flt(filename)
            % SEG2FLT converts a list of segments from a file to a list
            % of patches
            %
            %   flt=seg2flt(filename)
            
            
            fid=fopen(filename);  % open and count the number of columns
            line=strtrim(fgetl(fid));             
            while (strcmp('#',line(1)))
                line=strtrim(fgetl(fid));  % line ignores all shell comments and check for the first line that contains data
            end
            fclose(fid);
            
            [~,Vpl,x1,x3,width,d,wo,alphaw]=...
                textread(filename,'%u %f %f %f %f %f %f %f',...
                'commentstyle','shell');
                      
            fm=[];
            for k=1:length(x1)
                % list of patches for current segment
                flt=fgeom.flt2flt([x1(k);x3(k)],width(k),d(k),wo(k),alphaw(k));                                                   
                % list of patches for all segments
                fm=[fm;[flt,Vpl(k).*ones(size(flt,1))]];
            end
                       
        end
    end % methods (Static)
    methods(Static)
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %   convert segment definition to fault patches    %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function [fm] = loadFlt(filename)
            % LOADFLT loads the file and creates flt based on  
            % availability of Vpl
            %
            %   flt=loadFlt(filename)
            %
            %
            % where flt is a list of [N,n] patch geometry parameters
            
            fid=fopen(filename);  % open and count the number of columns
            line=strtrim(fgetl(fid));             
            while (strcmp('#',line(1)))
                line=strtrim(fgetl(fid));  % Open the file, get the first line
            end
            fclose(fid);
            
            [~,x1,x3,width,d]=...
                textread(filename,'%u %f %f %f %f',...
                'commentstyle','shell');
            fm = [x1,x3,width,d];
                                   
        end
    end % methods (Static)
end