# encoding: utf-8
'''
A tool for automatic transfer of minimal relevant data from the batch server to
local storage.

Created on 15 Jan 2019

@author: Eugen Rožić
'''
import os
import shutil
import time
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('server', type=str,
                    help='myriad or kathleen or lemon')
parser.add_argument('-u', '--user', type=str,
                    help='user name for the given server')
parser.add_argument('-s', '--search', type=str, default=r'growth',
                    help='the root search dir')
parser.add_argument('-f', '--full', action='store_true',
                    help='transfers the .dump files too')
parser.add_argument('--no_scp', action='store_true',
                    help="uses only the OS's sftp connection for transfer")
args = parser.parse_args()

if args.server == 'myriad':
    host = 'myriad.rc.ucl.ac.uk'
    user = 'ucapero' if args.user == None else args.user
    remote_local = r'/run/user/1000/gvfs/sftp:host={:s},user={:s}/'.format(host, user)
    remote_homedir = remote_local + 'lustre/home/' + user
    remote_rel_basedir = 'simulations/amyloids/data' #'Scratch/amyloids'
elif args.server == 'kathleen':
    host = 'kathleen.rc.ucl.ac.uk'
    user = 'ucapero' if args.user == None else args.user
    remote_local = r'/run/user/1000/gvfs/sftp:host={:s},user={:s}/'.format(host, user)
    remote_homedir = remote_local + 'lustre/home/' + user
    remote_rel_basedir = 'simulations/amyloids/data' #'Scratch/amyloids'
elif args.server == 'lemon':
    host = 'lemon0.biop.phys.ucl.ac.uk'
    user = 'erozic' if args.user == None else args.user
    remote_local = r'/run/user/1000/gvfs/sftp:host={:s},user={:s}/'.format(host, user)
    remote_homedir = remote_local + 'storage/users/' + user
    remote_rel_basedir = 'simulations/amyloids/data'
else:
    raise Exception('Unsupported server ({:s})!'.format(args.server))

remote_basedir = os.path.join(remote_homedir, remote_rel_basedir)
local_basedir = r'/media/data_ntfs/PhD/projects/amyloid MD simulations/data'

def remote_to_local(path):
    return os.path.join(local_basedir, os.path.relpath(path, remote_basedir))

def remote_to_scp(path):
    return r'{:s}@{:s}:"{:s}"'.format(user, host, os.path.relpath(path, remote_homedir))

def dq_esc(a_str):
    return '"' + a_str + '"'
    
def scp_transfer(remote_src, local_dest, attempts = 1):
    attempt = 1
    while attempt <= attempts:
        command = r'scp {:s} {:s}'.format(
            dq_esc(remote_to_scp(remote_src)), dq_esc(local_dest))
        exit_status = os.system(command)
        if exit_status == 0:
            break
        else:
            attempt += 1
    return exit_status

def sftp_transfer(remote_src, local_dest):
    print("Transfering {:s}".format(os.path.relpath(remote_src, remote_basedir)))
    shutil.copy(remote_src, local_dest)
    
def file_transfer(remote_src, local_dest, scp_attempts=1):
    if args.no_scp or scp_transfer(remote_src, local_dest, scp_attempts) != 0:
        sftp_transfer(remote_src, local_dest)

rootsearchdir = args.search
rootnode = os.path.join(remote_basedir, rootsearchdir)

if not os.path.exists(remote_basedir):
    raise Exception("Missing SFTP connection with host ({:s})!".format(host))
elif not os.path.exists(rootnode):
    raise Exception("The search root dir ({:s}) doesn't exist! Change the '-s' option...".format(
        os.path.relpath(rootnode, remote_homedir)))

start = time.time()
for node, dirs, files in os.walk(rootnode, topdown=True):
    node_name = os.path.basename(node)
    if node_name == 'old' or node_name+'.cfg' in files:
        # skip this dir and all subdirs (because either old or processed)
        dirs[:] = [] #this way because it has to remain the same object
        continue
    for filename in files:
        basename, ext = os.path.splitext(filename)
        if not ext == '.cfg':
            continue
        if basename in dirs: #this means that there exist some results...
            cfg_path = os.path.join(node, filename)
            dirs.remove(basename)
            srcdir = os.path.join(node, basename)
            print("Processing", os.path.relpath(srcdir, remote_basedir), "...")
            destdir = remote_to_local(srcdir)
            if not os.path.exists(destdir):
                os.makedirs(destdir)
            
            for result_filename in os.listdir(srcdir):
                if (result_filename.endswith('_last_dump') or
                    result_filename.endswith('_data') or
                    result_filename.endswith('_adsorbed') or
                    result_filename.endswith('.msd') or
                    result_filename.endswith('.lammps') or
                    (args.full and result_filename.endswith('.dump'))):
                    
                    result_filepath = os.path.join(srcdir, result_filename)
                    if not os.path.exists(remote_to_local(result_filepath)):
                        file_transfer(result_filepath, destdir)
                
                #TODO possibly delete some stuff (*.mol, *.dat, *.out, *.log)
            
            if not os.path.exists(os.path.join(destdir, filename)):
                file_transfer(cfg_path, destdir)
            
            if not os.path.exists(os.path.join(srcdir, filename)):
                shutil.move(cfg_path, srcdir)
            else:
                os.remove(cfg_path)
            
            print("... done!")
     
secs = time.time() - start
hrs = int(secs/3600); secs -= hrs*3600
mins = int(secs/60); secs -= mins*60
print("time: {:02d}:{:02d}:{:05.2f}".format(hrs, mins, secs))
