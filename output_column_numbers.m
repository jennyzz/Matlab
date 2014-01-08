function output_column_numbers(matrix,threshold)

lastLine = size(matrix,1);

for i=1:size(matrix,2)
    
    if (matrix(lastLine,i)>threshold)
        
        disp(sprintf('Column #%d exceeds thresold',i));
    end
   

end
