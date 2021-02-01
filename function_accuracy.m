function [accuracy] = function_accuracy(A,B)
%accuracy is a measure of mismatch between A and B. accuracy = 1 when A and
%B are linearly dependent. e.g. two identical images of different
%brightness.


prod = A.*B;
accuracy = sum(prod(:))/sqrt(sum(A(:).^2).*sum(B(:).^2));
end

