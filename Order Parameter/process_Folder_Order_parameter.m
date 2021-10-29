%Compute Order Parameter and generate results for each tif in a folder

%Get folder and parse input
folder = uigetdir;

path = dir(fullfile(folder, '*.tif'));
%Remove . and .. entries
path=path(~ismember({path.name},{'.','..'}));

results = cell(length(path),1);

outputs = cell(numel(results)+1,7);
header = ["File_Name","area", "aspect_ration", "avg_coherence", "Order_Parameter", ...
    "orientation", "Radial_Order"];
for k = 1:numel(header)
   outputs{1,k} = header(k); 
end

for i = 1:numel(results)
   file = fullfile(path(i).folder, path(i).name);
   im = import_tif(file);
   out = calculate_order_parameter(im);
   results{i} = out;
   outputs{i+1,1} = path(i).name;
   outputs{i+1,2} = out.area;
   outputs{i+1,3} = out.aspect_ratio;
   outputs{i+1,4} = out.avg_coherence;
   outputs{i+1,5} = out.S;
   outputs{i+1,6} = out.orientation;
   outputs{i+1,7} = out.S_radial;
end


% Convert cell to a table and use first row as variable names
T = cell2table(outputs(2:end,:),'VariableNames',header);
 
% Write the table to a CSV file
writetable(T,fullfile(folder,'Results.csv'))
