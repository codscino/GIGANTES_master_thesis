% 1. Create a sample workspace for demonstration
mySingleValue = 10;
myVector = [1, 2, 3, 4, 5];
myMatrix = rand(3,4);
myTable = table(rand(4,1), rand(4,1), 'VariableNames', {'Data1', 'Data2'});
mySimpleStruct = struct('name', 'John Doe', 'age', 30, 'isStudent', false);
myStructArray = struct('name', {'Jane Doe'; 'Peter Pan'}, 'age', {25; 12});
myNestedStruct = struct('Experiment', 'X123', 'data', mySimpleStruct);
myMixedStruct = struct('a', 1, 'b', {{'c', 'd'}}, 'c', mySimpleStruct);

% Save the workspace to a .mat file
save('myWorkspace.mat');
data = load('myWorkspace.mat');

% 2. Get fieldnames (each variable) in that struct
vars = fieldnames(data);

% 3. Open an output file
fid = fopen('workspace.md','w');

for k = 1:numel(vars)
    name = vars{k};
    value = data.(name);
    
    if isnumeric(value) && ismatrix(value)
        % Convert numeric array to a table for consistent processing
        T = array2table(value);
        T.Properties.VariableNames = matlab.lang.makeValidName(...
            compose("%s_%d", name, 1:width(T)) );
        
        fprintf(fid, "### %s\n\n", name);
        
        % Write table headers
        fprintf(fid, "|");
        for c = 1:width(T)
            fprintf(fid, " %s |", T.Properties.VariableNames{c});
        end
        fprintf(fid, "\n");
        
        % Write separator
        fprintf(fid, "|");
        for c = 1:width(T)
            fprintf(fid, " --- |");
        end
        fprintf(fid, "\n");
        
        % Write data rows
        for r = 1:height(T)
            fprintf(fid, "|");
            for c = 1:width(T)
                val = T{r,c};
                if isnumeric(val)
                    fprintf(fid," %g |", val);
                else
                    fprintf(fid," %s |", string(val));
                end
            end
            fprintf(fid, "\n");
        end
        
        fprintf(fid, "\n\n");
    elseif istable(value)
        % If it’s already a table, write it out
        T = value;
        fprintf(fid, "### %s\n\n", name);
        fprintf(fid, "|");
        for c = 1:width(T)
            fprintf(fid, " %s |", T.Properties.VariableNames{c});
        end
        fprintf(fid, "\n|");
        % Use single quotes for a char array and [] for concatenation
        fprintf(fid, [repmat(' --- |', 1, width(T))]);
        fprintf(fid, "\n");
        for r = 1:height(T)
            fprintf(fid,"|");
            for c = 1:width(T)
                cellVal = T{r,c};
                fprintf(fid," %s |", string(cellVal));
            end
            fprintf(fid,"\n");
        end
        fprintf(fid,"\n\n");
    elseif isstruct(value)
        % Handle struct data
        fprintf(fid, "### %s (struct)\n\n", name);
        
        struct_fields = fieldnames(value);
        
        % Write column headers from struct fields
        fprintf(fid, "|");
        for c = 1:numel(struct_fields)
            fprintf(fid, " %s |", struct_fields{c});
        end
        fprintf(fid, "\n");
        
        % Write separator row
        fprintf(fid, "|");
        for c = 1:numel(struct_fields)
            fprintf(fid, " --- |");
        end
        fprintf(fid, "\n");
        
        % Write data for each element in the struct array
        for r = 1:numel(value)
            fprintf(fid, "|");
            for c = 1:numel(struct_fields)
                field_name = struct_fields{c};
                val = value(r).(field_name);
                
                if isnumeric(val) || islogical(val)
                    fprintf(fid, " %g |", val);
                elseif ischar(val)
                    fprintf(fid, " %s |", val);
                elseif iscell(val)
                    fprintf(fid, " {[...]} |");
                elseif isstruct(val)
                    fprintf(fid, " [struct] |");
                else
                    fprintf(fid, " _%s_ |", class(val));
                end
            end
            fprintf(fid, "\n");
        end
        fprintf(fid, "\n\n");
    else
        % Skip other data types like cell arrays for now
        fprintf(fid, "### %s\n\n", name);
        fprintf(fid, "_(%s data not exported)_\n\n", class(value));
    end
end

fclose(fid);
disp('Workspace data has been exported to workspace.md');