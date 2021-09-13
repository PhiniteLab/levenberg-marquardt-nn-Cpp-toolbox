%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dynamic Neural Network Application with LM application
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Writing the dataset on the graph
t = [0:0.1:20];                           % input values
Y_general = sin(t);                       % target values

training_number = length(Y_general(1,:)); % number of training dataset

%% Input values are arranged.
X0 = ones(1,training_number);             % for bias input
X1 = t(1,:);                              % for first input

X_general = [X0;X1];                      % collective dataset

%% Writing the dataset on the graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% defining error parameters

errorNow = ones(training_number,1);  % error vector for now value
errorPre = ones(training_number,1);  % error vector for pre value

errorNowValue = sum(errorNow);       % error value for now value
errorPreValue = sum(errorPre);       % error value for pre value

epsilon = 0.000001;

%% defining error parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Neural Network parameters

number_of_input_layer_node = 2;
number_of_hidden_layer_node = 6;
number_of_output_layer_node = 1;

I = number_of_input_layer_node;
H = number_of_hidden_layer_node;
K = number_of_output_layer_node;

%% training parameters

iteration_max = 5000;        % for maximum iteration number
iteration = 0;               % for internal iteration number
detInvMatrices = 0;          % Jacobian determinant parameters
mu = 0.08;                   % learning parameter mu value
nnTrainingCondition = 1;     % neural network condition to be trained process


%% Neural Network parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% creating neural network structure

[W,W_previous,V,V_previous] = initialize_neural_network(H, K, I);

[z,y] = creating_activation_function(H, K, training_number);

%% creating neural network structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% general information related to the process

disp('Neural Network is started!');
disp('Basic information can be given by...');

disp('  ')
displayMessage = ['Input Layer Node: ',num2str(I), ' Hidden Layer Node :', num2str(H),...
    ' Output Layer Node: ',num2str(K)];
disp(displayMessage)

disp('  ')
disp('Training neural network is started in five seconds')

disp('  ')
pause(5);

nnOutputFileId = fopen('nnTrainInfo.txt','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% basic info transfer to file

fprintf(nnOutputFileId, strcat('Neural Network is started!','\n'));
fprintf(nnOutputFileId, strcat('Basic information can be given by...','\n'));
fprintf(nnOutputFileId, strcat(displayMessage,'\n'));
fprintf(nnOutputFileId, strcat('Training neural network is started in five seconds','\n'));

%% basic info transfer to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% general information related to the process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TRAINING PROCESS

%% comparing with the total error
while (nnTrainingCondition ~= 0)

    iteration = iteration + 1;
    errorPre = errorNow;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% training session is started
    JacobianTotal = [];

    for i = 1 : 1 : training_number

        % z value calculation
        z(:,i) = act_func_calc(X_general(:,i),W,H);

        % y value calculation
        y(:,i) = output_func_calc(z(:,i),V,K);

        JacobianTerm = phiLm(y(:,i), Y_general(:,i), z(:,i),X_general(:,i),V,K,H,I);
        
        JacobianTotal = [JacobianTotal;JacobianTerm];

        errorNow(i,1) = cost_function(y(:,i),Y_general(:,i));

    end

    %% training session is started
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    detInvMatrices = det(inv(JacobianTotal'*JacobianTotal + mu*eye(H*I + K*(H+1))));
 
    coeffUpdate = (JacobianTotal'*JacobianTotal + mu*eye(H*I + K*(H+1)))\(JacobianTotal'*errorNow);
    
    [W_new,V_new] = phiLmUpdate(W,V,coeffUpdate,I,H,K);

     W = W_new;
     V = V_new;

    errorNowValue = sum(errorNow)/training_number;
    errorPreValue = sum(errorPre)/training_number;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (errorPreValue - errorNowValue) > 0

        internalAssessmentMu = (errorPreValue - errorNowValue);

        if abs(internalAssessmentMu) > (1e-2/training_number)

            mu = mu + mu*0.01;

        else

            mu = mu - mu*0.01;

        end

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% displaying the whole results

    if (mod(iteration,10) == 0)
    
        displayMessage = ['Error: ',num2str(errorNowValue),' Iteration: ',...
            num2str(iteration), ' Jacobian Check: ',num2str(detInvMatrices),' Mu: ',num2str(mu)];

        fprintf(nnOutputFileId,strcat(displayMessage,'\n'));
        
        disp(displayMessage)

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% nn condition update
    
    nnTrainingCondition = (errorNowValue > epsilon) && (iteration < iteration_max)...
    && (mu > 1e-7) && (detInvMatrices < 1e200);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

end

%% TRAINING PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST PROCESS

z_test = zeros(size(z));
y_test = zeros(size(y));

for i = 1 : 1 : training_number

    % z value calculation
    z_test(:,i) = act_func_calc(X_general(:,i),W,H);

    % y value calculation

    y_test(:,i) = output_func_calc(z(:,i),V,K);
    
end

figure
Y_plot = Y_general;

plot(Y_plot)

y_model_plot = y_test;

hold on
plot(y_model_plot)

%% TEST PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(nnOutputFileId)


%% Dynamic Neural Network Application
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%