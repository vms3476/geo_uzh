#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:14:50 2010

@author: Rémi Flamary
"""


import glob
import getopt,sys,os,psutil,argparse
from time import strftime,strptime
from datetime import datetime,timedelta


dir_proc=os.path.expanduser("~/.remtime/procs/")

conf_dir=os.path.expanduser("~/.remtime/")

class TaskInfo:
    """ Class that handles each task"""
    task_name=""
    time_start=datetime.now()
    percent_done=0;
    time_percent_done=datetime;
    delta_till_next=datetime;
    est_finite_time=datetime;
    est_remaining_time=timedelta;
    last_comment=""
    tot_time=0;
    nbcall=0
    pid=0

    def __init__(self,taskname):
        self.task_name=taskname
        if len(taskname):
            self.load_task(taskname)
    
    def get_rem_time(self):

        if self.percent_done>0 and self.percent_done<1:
    

            # compute predicted finish time
            frac = 1/self.percent_done
            dif_time=self.time_percent_done-self.time_start
            #print "diff",dif_time
            #tot_time=frac*dif_time.total_seconds();
            self.tot_time=int((dif_time.microseconds + (dif_time.seconds + dif_time.days * 24 * 3600) * 10**6) / 10**6)

            tot_time=frac*self.tot_time
            #print "tot_time",tot_time
            tot_time=int(tot_time)

            tot_time=timedelta(seconds=tot_time);
            #print "tot_time",tot_time
            self.est_finite_time=self.time_start+tot_time;
            if datetime.now()>self.est_finite_time:
                self.est_remaining_time= timedelta(seconds=0)
            else:
                self.est_remaining_time=self.est_finite_time-datetime.now();

            self.est_remaining_time=timedelta(self.est_remaining_time.days,self.est_remaining_time.seconds);

            # compute the delta before the next one
            dif_time=self.time_percent_done-self.time_start

            tot_time=self.tot_time/self.nbcall;
            tot_time=int(tot_time)
            tot_time=timedelta(seconds=tot_time);
            self.delta_till_next=tot_time

        return self.est_remaining_time

    def get_time_warning(self):
        if datetime.now()>(self.time_percent_done+self.delta_till_next):
            return "[!]"
        else:
            return ""
    

    def load_task(self,task):
        self.task_name=task
        f = open(dir_proc+task, 'r')
        lines=f.readlines()
        l1=lines[0]
        lf=lines[-1]
        l1s=l1.split("\t")
        lfs=lf.split("\t")
        self.nbcall=len(lines)-1;
        
        self.percent_done=float(lfs[0])
        self.pid=int(lfs[1]);
        self.time_start=datetime.strptime(l1s[2],"%Y-%m-%d %H:%M:%S")
        #print "Starting time: ",self.time_start.strftime("%Y-%m-%d %H:%M:%S")
        self.time_percent_done=datetime.strptime(lfs[2],"%Y-%m-%d %H:%M:%S")
        self.last_comment=lfs[3]
        #print "Last percent time: ",self.time_percent_done.strftime("%Y-%m-%d %H:%M:%S")
        self.get_rem_time()
        #else: 


        

    def get_pid(self):
        return self.pid

    #def percent_done(self):
    #    return self.percent_done
    


    def print_task(self):
        print "Task Name: ",self.task_name
        print "Percent done: ",self.percent_done*100
        print "Starting time: ",self.time_start.strftime("%Y-%m-%d %H:%M:%S")
        if self.percent_done>0 and self.percent_done<1:
           print "Est. remaining time: ",self.est_remaining_time,self.get_time_warning()
           print "Est finish time: ",self.est_finite_time.strftime("%Y-%m-%d %H:%M:%S"),self.get_time_warning()
        elif self.percent_done==0 :
            print "No prevision yet"
        elif self.percent_done==1:
            print "Task Completed!"
        if len(self.last_comment)>0:
            print "Last Comment: ",self.last_comment

        # process info
        if self.percent_done<1:
            p=psutil.Process(self.pid)
            print "Cmdline:",p.cmdline
            print "CPU times:",p.get_cpu_times()
            print "CPU %:",p.get_cpu_percent()
        print

    def get_task_info_txt(self):
        res= 'Name:\n' + self.task_name+ "\n"
        res+= "Percent done: " + str(self.percent_done*100) + "\n"
        res+= "Starting time: \n" + self.time_start.strftime("%Y-%m-%d %H:%M:%S")+ "\n"
        if self.percent_done>0 and self.percent_done<1:
           remtime = str(self.est_remaining_time) + self.get_time_warning()       
           finishtime=self.est_finite_time.strftime("%Y-%m-%d %H:%M:%S") + self.get_time_warning()
           res+= "Est. remaining time: \n" + remtime + "\n"
           res+= "Est finish time: \n" + finishtime + "\n"
        elif self.percent_done==0 :
            res+= "No prevision yet"
        elif self.percent_done==1:
            res+= "Task Completed!"   
        return res
        

#        print "Task Name: ",self.task_name
#        print "Percent done: ",self.percent_done*100
#        print "Starting time: ",self.time_start.strftime("%Y-%m-%d %H:%M:%S")

#        if len(self.last_comment)>0:
#            print "Last Comment: ",self.last_comment
#
#        # process info
#        if self.percent_done<1:
#            p=psutil.Process(self.pid)
#            print "Cmdline:",p.cmdline
#            print "CPU times:",p.get_cpu_times()
#            print "CPU %:",p.get_cpu_percent()
#        print
        
    def print_short(self):
        if self.percent_done>0 and self.percent_done<1:
            print self.percent_done*100,"%\t",self.task_name,"\t",self.get_rem_time(),self.get_time_warning()
        elif self.percent_done==0 :
            print self.percent_done*100,"%\t",self.task_name,"\t","Not Available"
        elif self.percent_done==1:
            print self.percent_done*100,"%\t",self.task_name,"\t","Finished"

def get_str_pid(m):
    """Get the PPID of PPPID, used to store the task PID"""
    if m:
        p=psutil.Process(os.getppid())
        #print p
        str_pid=str(p.ppid)
    else:
        str_pid=str(os.getppid())
    return str_pid

def new(task,message,matlab):
    """Create a new task with agivem message (get the PPPID if matlab)"""
    #print "Adding task"
    #print task
    if len(task)>0:
        ctime=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f = open(dir_proc+task, 'w')
        str_pid=get_str_pid(matlab)
        f.write("0\t" + str_pid + "\t" + ctime + "\t"+ message+ "\n" )
        f.close()

def update(task,done,message,matlab):
    """Update the given task with the value done (0,1)""" 
    #print "Updating task"    
    #print "Adding task"
    #print task
    if len(task)>0:
        ctime=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f = open(dir_proc+task, 'a')
        str_pid=get_str_pid(matlab)
        f.write(str(done) + "\t"+ str_pid + "\t" + ctime + "\t"+ message+ "\n" )
        f.close()

def kill(task):
    """Kill the given task (from task list and process)"""
    print "Killing task:",task
    if len(task)>0:
        # killing task if still running
        if os.path.isfile(dir_proc+task):
            info=get_task_infos(task)
            os.remove(dir_proc+task)
            try:
                p=psutil.Process(info.get_pid())
                if p.is_running():
                    p.kill()
                # deleting file
            except psutil.error.NoSuchProcess:
                print "System Process does not exists"
        else:
            print "Task does not exists"

def killexp(exp,force=0):
    """ Kill all the tasks whose name correspond to the given regexp"""
    print "Killing tasks with expression:",exp
    none=False
    if len(exp)>0:
        for task in glob.glob(dir_proc+exp):
            if force:
                task=os.path.basename(task)
                kill(task)
            elif not none:                                 
                print "Delete task : ",os.path.basename(task),"? [y/n/all/none]",
                answer = raw_input()
                if answer == 'y':
                    task=os.path.basename(task)
                    kill(task)
                elif answer == 'all':
                    force=True
                elif answer == 'none':
                    none=True
    

def cleantask(task):
    """Clean the given task from the task list"""
    if len(task)>0:
        # killing task if still running
        if os.path.isfile(dir_proc+task):
            info=get_task_infos(task)
            os.remove(dir_proc+task)
        else:
            print "Task does not exists"
            
def get_task_list():
    res=os.listdir(dir_proc)
    res.sort(key=lambda x: x.lower())
    return  res

def get_percent_task_list():
    res=os.listdir(dir_proc)
    res.sort(key=lambda x: x.lower())
    for i in range(len(res)):
        info=get_task_infos(res[i])
        res[i]= "{0:2.1f}".format(info.percent_done*100) + "%  " + res[i]
    return res         
            
def clean():
    """Clean all the finished tasks from the task list"""
    #print "Killing finisehd task"
    listtasks=os.listdir(dir_proc)
    for task in listtasks:
        info=TaskInfo(task)
        #print info.get_txt()
        if info.percent_done==1:
            print "Cleaning task:",task
            cleantask(task)

def view(task):
    """Print informations about the given task"""
    info=get_task_infos(task)
    #print info.get_txt()
    info.print_task();

def vlist():
    """Print a list of tasks"""
    print "Done \tTask \tRem Time" 

    listtasks=os.listdir(dir_proc)
    listtasks.sort(key=lambda x: x.lower())
    for task in listtasks:
        info=TaskInfo(task)
        #print info.get_txt()
        info.print_short();

def get_task_infos(task):
    # get task infios     
    info=TaskInfo(task)
    return info
    
def check_install():
    
    if not os.path.exists(conf_dir):
        print("Creating config files")
        os.mkdir(conf_dir)  
        os.mkdir(dir_proc)  

def main(argv):  

 # check if installed properly
    check_install()
    
    
    parser = argparse.ArgumentParser(prog='remtime',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''remtime is a tool for tracking advancements of long tasks. You can list them, estimate the remaining time for each task and stop them from command line''',
    epilog=''' 
 Examples:
    Creating a task name "test_task" with message "Intitialization":
       remtime new "test_task" -m "Intitialization"

    Updating the task "test_task" after 10% done:
       remtime update "test_task" 0.10

    Visualizing the task list and the percent done:
       remtime list 

    Visualizing task "test_task" and see estimated finish time:
       remtime view "test_task"

    Cleaning the tasks (delete finished ones):
       remtime clean 

    Deleting the task "test_task" (kills the process too):
       remtime kill "test_task" 

    
For any comments contact the author''')   
    
    #parser.add_argument('task', metavar='task', type=str, nargs=1,
    #               help='the task that should be performed',action="store")

               
    #parser.add_argument('-b','--bibfile',dest='bibfile',type=str, help="change from default bib file",default='')
    #parser.add_argument('-c','--congigfile',dest='configfile',type=str, help="change from default config file",default='')
                           
                   
    subparsers = parser.add_subparsers(help='task that will be performed',dest="task")
    
    parser_new = subparsers.add_parser('new', help='create a new task ')
    parser_new.add_argument('taskname', type=str, help='unique task name (string)')
    parser_new.add_argument('-m','--message', type=str, help='optionnal message stored for the task',default='')
    parser_new.add_argument('--matlab', action='store_true', help='use if calling from matlab (selects the father of the callinc process)')  

    parser_update = subparsers.add_parser('update', help='update the task state (percent done)')
    parser_update.add_argument('taskname', type=str, help='unique task name (string)')
    parser_update.add_argument('done', type=float, help='state of the task ( between 0.0-1.0)')
    parser_update.add_argument('-m','--message', type=str, help='optionnal message stored for the task',default='')
    parser_update.add_argument('--matlab', action='store_true', help='use if calling from matlab (selects the father of the callinc process)')  


    parser_view = subparsers.add_parser('view', help='view information about task')
    parser_view.add_argument('taskname', type=str, help='unique task name (string)')
    
    parser_list = subparsers.add_parser('list', help='list the current tasks')
    
    parser_clean = subparsers.add_parser('clean', help='clean the list of tasks (remove finished)')

    parser_kill = subparsers.add_parser('kill', help='kill the selected task')    
    parser_kill.add_argument('taskname', type=str, help='unique task name or expression for several tasks')    
    parser_kill.add_argument('-f','--force', action='store_true', help='kill tasks silently')    
    
    args= parser.parse_args()   
               
    todo= args.task                     
    if todo == "new":
        new(args.taskname,args.message,args.matlab)
    elif todo == "update":
        update(args.taskname,args.done,args.message,args.matlab)
    elif todo == "kill":
        killexp(args.taskname)
    elif todo == "view":
        view(args.taskname) 
    elif todo == "list":
        vlist() 
    elif todo == "clean":
        clean()
        
        
if __name__ == "__main__":
    main(sys.argv[1:])
	
