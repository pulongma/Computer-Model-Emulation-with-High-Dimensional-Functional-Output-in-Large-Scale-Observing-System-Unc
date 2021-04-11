function [output_sorted, input_sorted, ind] = sort_data(output, input, sort_dim)

n = size(output, 1);
q = size(output, 2);

[output_new, ind] = sort(input(:,sort_dim));

output_sorted = output(ind, :);
input_sorted = input(ind, :);



end 

