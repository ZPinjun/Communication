function data = func_PerSheetRead(filename)
% Obtain sheet names
sheets = sheetnames(filename);
% Create a cell to store all sheet contents
data = cell(length(sheets), 1);
% Iterate sheets
for i = 1:length(sheets)
    data{i} = readmatrix(filename, 'Sheet', sheets{i});
end
end
