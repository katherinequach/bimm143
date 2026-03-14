## Basic Unixs

Some important file system commands include

pwd: Print working directory
ls: List files and folders
sd: Change directory
mkdir: Make a new directory
rm: Delete files and directories **Careful**
nano: A very basic text editor that is always available
less: To view/read text files page by page (pager program)


My AWS instance:
ssh -i ~/Downloads/bimm143_kq.pem ubuntu@ec2-44-255-46-28.us-west-2.compute.amazonaws.com

To copy from my AWS instance
scp -i ~/Downloads/bimm143_kq.pem ubuntu@ec2-44-255-46-28.us-west-2.compute.amazonaws.com:~/work/results.tsv .


## Class 17 AWS Instance Address

ssh -i "bimm143_kq.pem" ubuntu@ec2-54-244-192-79.us-west-2.compute.amazonaws.com