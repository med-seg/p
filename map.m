function [F] = map(X, Y)

% This function extract parameter vales corresponding to rows X from table Y

X = abs(round(X));
Y = abs(round(Y));

if(X <= 0)
X = 1;
end

if(Y <= 0)
Y = 1;
end
   
P = X >= length(Y)
X(P) = 0;
        
F = Y(X+1)
F(P) = 0;

end

