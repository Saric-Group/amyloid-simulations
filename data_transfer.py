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

parser.add_argument('user_host', type=str,
                    help='user@host for SCP')
parser.add_argument('-s', '--search', type=str, default=r'growth',
                    help='the root search dir')
args = parser.parse_args()

scp_prefix = args.user_host
user, host = scp_prefix.split('@')
remote_local = r'/run/user/1000/gvfs/sftp:host={:s},user={:s}/'.format(host, user)
if host.startswith('myriad'):
    remote_homedir = remote_local + 'lustre/home/ucapero'
    remote_rel_basedir = r'Scratch/amyloids'
elif host.startswith('lemon'):
    remote_homedir = remote_local + 'storage/users/erozic'
    remote_rel_basedir = r'simulations/amyloids/data'
else:
    raise Exception('Unsupported host ({:s})!'.format(host))

remote_basedir = os.path.join(remote_homedir, remote_rel_basedir)
local_basedir = r'/media/data_ntfs/PhD/simulation stuff/amyloids/data'

def remote_to_local(path):
    return os.path.join(local_basedir, os.path.relpath(path, remote_basedir))

def remote_to_scp(path):
    return r"{:s}:'{:s}'".format(scp_prefix, os.path.relpath(path, remote_homedir))

def qesc(str):
    if not '"' in str:
        return r'"{:s}"'.format(str)
    elif not "'" in str:
        return r"'{:s}'".format(str)
    else:
        raise Exception("Can't escape a string that contains both a single and double quotation mark!")
    
def scp_transfer(remote_src, local_dest, attempts = 3):
    attempt = 1
    while attempt <= attempts:
        exit_status = os.system(r'scp {:s} {:s}'.format(
            qesc(remote_to_scp(remote_src)), qesc(local_dest)))
        if exit_status == 0:
            break
        else:
            attempt += 1
    return exit_status

rootsearchdir = args.search
rootnode = os.path.join(remote_basedir, rootsearchdir)

start = time.time()
for node, dirs, files in os.walk(rootnode):
    if os.path.basename(node) == 'old':
        del dirs
        continue
    for filename in files:
        basename, ext = os.path.splitext(filename)
        if not ext == '.cfg':
            continue
        if os.path.basename(node) == basename:
            # a dir containing a .cfg with the same name is a cue for it having been
            # processed and to be skipped
            #del dirs #(shouldn't contain nested dirs, only files)
            break
        if basename in dirs: #this means that there exist some results...
            cfg_path = os.path.join(node, filename)
            dirs.remove(basename)
            srcdir = os.path.join(node, basename)
            print "Processing", os.path.relpath(srcdir, remote_basedir), "..."
            destdir = remote_to_local(srcdir)
            if not os.path.exists(destdir):
                os.makedirs(destdir)
            
            error = 0
            for result_filename in os.listdir(srcdir):
                if (result_filename.endswith('_last_dump') or
                    result_filename.endswith('_cluster_data')):
                    
                    result_filepath = os.path.join(srcdir, result_filename)
                    if os.path.exists(remote_to_local(result_filepath)):
                        continue
                    error += scp_transfer(result_filepath, destdir)
                
                #TODO possibly delete some stuff (*.mol, *.dat, *.out, *.log)
            
            if not os.path.exists(os.path.join(destdir, filename)):
                error += scp_transfer(cfg_path, destdir)
            
            if error == 0: #if everything was successful
                shutil.move(cfg_path, srcdir)
                print "... done!"
            else:
                print "... not fully processed! Check manually."
     
secs = time.time() - start
hrs = int(secs/3600); secs -= hrs*3600
mins = int(secs/60); secs -= mins*60
print "time: {:02d}:{:02d}:{:05.2f}".format(hrs, mins, secs)
