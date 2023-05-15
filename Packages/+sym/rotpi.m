function [state] = rotpi(domain,state)

s3 = size(state,3);
if(s3 == 1)
    state = flip(flip(state,1),2);
elseif(s3 == 2)
    state(:,:,1) = -flip(flip(state(:,:,1),1),2);
    state(:,:,2) = -flip(flip(state(:,:,2),1),2);
end

end

