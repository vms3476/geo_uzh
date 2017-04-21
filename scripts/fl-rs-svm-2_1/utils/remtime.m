function remtime(task,done,message)
% matlab function for the use of remtime

if nargin<3
    message='';
end

if nargin<2
    done=0;
end


if done==0
    [status, result] =system(['remtime new "' task '" -m "' message '" --matlab']);
else
    [status, result] =system(['remtime update "' task '" ' num2str(done) ' -m "' message '" --matlab']);
end
