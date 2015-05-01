classdef parobj < handle
    %PAROBJ Parameter object. Methods to set parameter values, types and
    % generate sets of parameters to sweep over using the method parsweep.
    %
    % See example of basic use in example_parobj.m.
    %
    % Jakob S. Joergensen (jakj@dtu.dk), 2014.
    
    
    properties (SetAccess = private)
        values = {};
        names  = {};
        types  = {};
        array  = [];
        stub   = 'sim';
        is_sync_values_array = true;
    end
    
    methods

        function setValues(PO,values)
            % For specifying the list of values for each parameter. The
            % input values must be a cell array with the length equal to
            % the number of parameters and each cell contain a vector with
            % the desired values for the given parameter.
            
            % If values is unset, set default Names and Types
            if isempty(PO.values)
                num_pars = length(values);
                PO.types = cell(1,num_pars);
                for k = 1:num_pars
                    PO.types{k} = '%f';
                end
            end
            PO.values = values;
            PO.is_sync_values_array = false;
        end
        
        
        function setNames(PO,names)
            % For specifying the names of each parameter. The
            % input names must be a cell array with the length equal to
            % the number of parameters and each cell contain a string with
            % the desired name for the given parameter.
            PO.names = names;
        end
        
        
        function setTypes(PO,types)
            % For specifying the datatype for each parameter. The
            % input types must be a cell array with the length equal to
            % the number of parameters and each cell contain a small string
            % specifying the datatype in fprintf format, eg '%d' for an
            % integer.
            PO.types = types;
        end
        
        function setStub(PO,stub)
            % For specifying the initial part of the automatically
            % generated filename. The input stub must be a string.
            PO.stub = stub;
        end
        
        
        function buildArray(PO)
            % Reads the specified values of the parameters and creates
            % array with all possible combinations in the rows.
            num_pars = length(PO.values);
            PC = cell(num_pars,1);
            if num_pars == 1
                PC{1} = PO.values{1}(:);
            else
                [PC{:}] = ndgrid(PO.values{end:-1:1});
            end
            PO.array = zeros(numel(PC{1}),num_pars);
            for k = 1:num_pars
                PO.array(:,k) = PC{num_pars - k + 1}(:);
            end
            PO.is_sync_values_array = true;
        end
        
        
        function setArray(PO,array)
            % Manually set the parameter object's array to the input array.
            % Note that this makes the parameter object's array
            % inconsistent with what might be specified in the values.
            PO.array = array;
            PO.is_sync_values_array = false;
        end
        
        
        function export(PO,filename)
            % Exports the parameter object to a text file with all
            % parameter sets listed.
            fid = fopen(filename,'w');
            
            num_pars = length(PO.values);
            
            % Write the first line with names of the parameters
            fprintf(fid,['#', repmat(' %s',1,num_pars),'\n'],PO.names{:});
            
            % Write the second line with types of the parameters
            fprintf(fid,['#', repmat(' %s',1,num_pars),'\n'],PO.types{:});
            
            % Make the format string for printing each parset to file.
            printstr = [sprintf('%s ',PO.types{:}),'\n'];
            
            % Write to the file
            fprintf(fid, printstr, PO.array');
            fclose(fid);
        end
        
        
        function import(PO,filename)
            % Read in a text file with parameters sets into a parobj.
            fid = fopen(filename,'r');
            
            C = textscan(fid,'# %s %s %s','CollectOutput',1);
            
            names = C{:}(1,:);
            types = C{:}(2,:);
            PO.setNames(names);
            PO.setTypes(types);
            num_pars = size(names,2);
            
            fmtstr = sprintf(repmat('%s ',1,num_pars),types{:});
            PO.array = fscanf(fid,fmtstr,[num_pars,inf])';
            
            PO.is_sync_values_array = false;
            
            
        end
        
        function fname = buildFormatString(parobj)
            % Build the automatic filename string from the parameter
            % values.
            
            fname = parobj.stub;
            
            num_pars = length(parobj.names);
            
            if isempty(parobj.types)
                parobj.types = cell(num_pars,1);
                for k = 1:num_pars
                    parobj.types{k} = '%f';
                end
            end
            
            
            for k = 1:num_pars
                fname = [fname, '_',parobj.names{k},'_',parobj.types{k}];
            end
        end
        
        function [] = parsweep(parobj, sim_func, pass_savefile)
            % Main method which sweeps over all rows of parameters present
            % in the object's array field, passing them as inputs to the
            % sim_func. The boolean input pass_savefile should be true if
            % the sim_func expects one additional input compared to the
            % number of parameters; this last input will be a string
            % generated by parsweep to be used for naming output files
            % systematically from each of the specific parameter sets that
            % the sim_func is called with during the sweep. If
            % pass_savefile is set to false, then no additional input is
            % expected and the sim_func should itself set up the names for
            % file(s) it saves.
            
            if nargin < 3
                pass_savefile = 'true';
            end
            
            [num_parsets, num_pars] = size(parobj.array);
            
            if isempty(parobj.names)
                parnames = cell(num_pars,1);
                for k = 1:num_pars
                    names{k} = sprintf('p%d',k);
                end
                parobj.setNames(names);
            end
            
            if isempty(parobj.types)
                types = cell(num_pars,1);
                for k = 1:num_pars
                    types{k} = '%e';
                end
                parobj.setTypes(types);
            end
            
            if  isempty(parobj.stub)
                parobj.setStub('sim');
            end
            
            % Save filename (not path): Only dependent on parameters
            savefileformatstr = parobj.buildFormatString();
            
            fprintf(['---------------------------------------------------------\n',...
                '---------------------------------------------------------\n',...
                '|   PARSWEEP BEGIN                                       \n',...
                '---------------------------------------------------------\n',...
                '---------------------------------------------------------\n']);
            % Main loop: Extract a row of parameters and pass to sim_func.
            for k = 1:num_parsets
                fprintf([...
                    '---------------------------------------------------------\n',...
                    '|   Running sim %d of %d...\n', ...
                    '---------------------------------------------------------\n'],...
                    k, num_parsets)
                pars = parobj.array(k,:);
                cellpars = num2cell(pars);
                savefile = sprintf(savefileformatstr, cellpars{:});
                if pass_savefile
                    feval(sim_func, ...
                        cellpars{:}, ...
                        savefile)
                else
                    feval(sim_func, ...
                        cellpars{:})
                end
                fprintf('DONE.\n')
            end
            fprintf(['---------------------------------------------------------\n',...
                '---------------------------------------------------------\n',...
                '|   PARSWEEP DONE                                        \n',...
                '---------------------------------------------------------\n',...
                '---------------------------------------------------------\n']);
        end
        
        function outcell = loadResults(parobj, ...
                loader, ...
                resfilepre, ...
                resfilepost, ...
                varargin)
            % Load results in files generated by parsweeping the parobj.
            %
            % Inputs:
            % - loader: Name of function to load (and perhaps process) single file.
            % - parobj: The parameter object desired to be loaded.
            % - resfilepre: String to prepend to the results file names,
            % typically this could be the path.
            % - resfilepost: String to append to the results file names.
            % - varargin: Any additional arguments are passed directly to loader
            
            % First build the filename format string to use in loading
            fileformatstr = [resfilepre, parobj.buildFormatString(),resfilepost];
            
            % Number of parameter sets and number of variables returned by loader
            num_parsets = size(parobj.array,1);
            num_loaded_variables = nargout(loader);
            
            % Set up cell array to hold the loaded results. Each row is for a single
            % parset/datefile, each column is for each output from the loader file.
            outcell = cell(num_parsets, num_loaded_variables);
            
            % Loop over the parsets: Load data from file given by the current parset
            % using loader.
            for k = 1:num_parsets
                fprintf('Loading %d of %d... ',k,num_parsets);
                pars = parobj.array(k,:);
                cellpars = num2cell(pars);
                savefile = sprintf(fileformatstr, cellpars{:});
                [outcell{k,:}] = feval(loader,savefile,varargin{:});
                fprintf('DONE\n');
            end
        end
        
        function cube = cubify(parobj,vec)
            % Turn a vector with an entry for each parset into a hyber-cube
            % with the number of dimensions equal to the number of
            % parameters and the length in each dimension equal to the
            % number of values of each parameter.
            
            dims = zeros(size(parobj.values));
            for k = 1:length(dims)
                dims(k) = length(parobj.values{k});
            end
            cube = reshape(vec,dims(end:-1:1));
            cube  = permute(cube,length(size(cube)):-1:1);
        end
        
    end
end
