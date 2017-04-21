function multi_remtime(task,depth,done,message)
% matlab function for the use of remtime
% for initialisation done=0, message=nbloop


if nargin<4
    message='';
end

global remtime_nbs
global remtime_done

if done==0
    
    remtime_nbs(depth)=message;
    remtime_done(depth)=0;
    
else
    remtime_done(depth)=done;
    done=compute_percent(1);
    remtime(task,done,message);
       
end


end

function done=compute_percent(p)
% recursive function for percent computing


global remtime_nbs
global remtime_done

if p>length(remtime_nbs)
   done=0;    
else
    
    done=remtime_done(p)/remtime_nbs(p)+1/remtime_nbs(p)*compute_percent(p+1);
    
    if (remtime_done(p)==remtime_nbs(p)) && p>1
        remtime_done(p)=0;
        remtime_done(p-1)=remtime_done(p-1)+1;
    end
    
    
end

end




