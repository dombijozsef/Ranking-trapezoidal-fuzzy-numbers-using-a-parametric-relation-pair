clear

%% FUZZY NUMBERS
% There are the parameters of the input trapezoidal fuzzy numbers 
fuzzy_sets = [
    33    65    68    74
    55    58    62    73
    67    76   125   130
    42    71    87   111
    72    85    88    91
    65    68    72    92
    57    72    77    85
    28    31    34    61
     2     5    20    66
    11    70   73    76];

% Here we get delta for wich the relation is transitive
delta = interval_preference_relation_class.get_intransitive_triples(fuzzy_sets)

% Here we rank the fuzzy numbers
[orderedFuzzySets,orderedFuzzyIds] = interval_preference_relation_class.order_fuzzy_numbers(delta,fuzzy_sets);

% Here we plot the original and the ordered numbers
interval_preference_relation_class.plot_fuzzy_numbers(fuzzy_sets,1:size(fuzzy_sets,1),'original_set.eps','original');
interval_preference_relation_class.plot_fuzzy_numbers(orderedFuzzySets,orderedFuzzyIds,'ordered_set.eps','ordered');
