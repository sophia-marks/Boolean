% "one_network_iteration_rules.m"
% Inputs:
% inputvector = the initial state, a vector
% parents = neighborhood of each node, obtained with a different function
% rule = rule array (result from binary_rule)

function [outputvector,outputlabels] = one_network_iteration_rules_StateLabels(inputvector,parents,rule)

for index=1:length(inputvector)
    t = inputvector(parents{1,index});
    % Convert the vector of parent values to a string
    T = num2str(t);
    % Remove the spaces from the string
    T = strrep(T,'1 ','1');
    T = strrep(T,'0 ','0');
    T = strrep(T,'1 ','1');
    T = strrep(T,'0 ','0');
    outputlabels(index) = deblank(string(T));
    % Find the correct result, using the decimal representation of the
    % parent values as an index for the rule array
    outputvector(index) = rule(2^length(t)-bin2dec(T)); %#ok<AGROW>
end
