% idinput() can also be used to generate m-sequences
function[pilot] = ruo_pilot_gen(pri_poly)
    n = length(pri_poly); % 移位寄存器长�?
    N = 2^n-1; % 伪随机码的周�?
    register = [zeros(1,n-1) 1];
    pilot = zeros(1,N);
    for i = 1:N
        newregister = mod(sum(pri_poly.*register),2);
        register(2:n) = register(1:(n-1));
        register(1) = newregister;
        pilot(i) = register(n);
    end
end

% function[pilot] = ruo_pilot_gen(pri_poly)
%     n = length(pri_poly); % 移位寄存器长�?
%     N = 2^n-1; % 伪随机码的周�?
%     register = [zeros(1,n-1) 1];
%     pilot = zeros(1,N);
%     pilot(1) = register(n);
%     for i = 2:N
%         newregister = mod(sum(pri_poly.*register),2);
%         register(2:n) = register(1:(n-1));
%         register(1) = newregister;
%         pilot(i) = register(n);
%     end
% end

