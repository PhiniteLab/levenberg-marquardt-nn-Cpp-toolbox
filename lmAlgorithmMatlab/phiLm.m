function [JacobianTerm] = phiLm(nnOutput,trainingOutput, activationOutput, inputData, V_weights, K, H, I)

J_int = [];

for k = 1 : 1 : K

    for h = 1 : 1 : H + 1

        J_int = [J_int,(trainingOutput(k,1) - nnOutput(k,1))*activationOutput(h,1)];

    end

end

for h = 1 : 1 : H

    for in = 1 : 1 : I

        J_int = [J_int,(trainingOutput(:,1) - nnOutput(:,1))'*V_weights(:,h)*activationOutput(h,1)*(1 - activationOutput(h,1))*inputData(in,1)];

    end

end

JacobianTerm = J_int;

end
