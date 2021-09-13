function [W,W_previous,V,V_previous] = initialize_neural_network(number_of_hidden_layer_node, number_of_output_layer_node, number_of_input_layer_node)

    I = number_of_input_layer_node;
    
    H = number_of_hidden_layer_node;
    
    K = number_of_output_layer_node;

    W = randi([-10 10],H,I)./100;

    W_previous = randi([-10 10],H,I)./100;

    V = randi([-10 10],K,H+1)./100;

    V_previous = randi([-10 10],K,H+1)./100;


end