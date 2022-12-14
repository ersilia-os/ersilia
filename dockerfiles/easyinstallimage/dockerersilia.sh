#Check if the ~/ersiliafiles exists. If not, create one.
echo "Welcome to ersilia docker image quickstart!"
echo "This script will fetch and run ersilia docker image."
echo "It will also create a directory called \"ersiliafiels\" in your home directory."
echo "This is the location where ersilia stores all the data."
#read -p "Do you wish to continue?"
#if [[ ! $REPLY =~ ^[Yy]$ ]]
#then
#        exit
#fi

cd ~ || exit
if [ ! -d $(pwd)/ersiliafiles ]
then
        echo "Directory ~/ersiliafiles doesn't exist. Creating it now."
        mkdir $(pwd)/ersiliafiles
fi
#Docker startup
if ! docker ps|grep -q "ersilia"
then
        sudo docker run --name ersilia -d --mount type=bind,source=$(pwd)/ersiliafiles,target=/root/eos ersiliaos/ersilia
else
        echo "Ersilia docker container already started"
fi
#Shell detection
if echo $SHELL|grep -q "zsh"; then
        rcfile=~/.zshrc
elif echo $SHELL|grep -q "bash"; then
        rcfile=~/.bashrc
elif ps -p $$|grep -q "sh"; then
        rcfile=~/.shrc
else    echo "Shell unknown. Please add the ersilia alias manually to your rc file."
        echo "alias ersilia="docker exec -it ersilia /bin/bash""
fi
echo  "Your rc file has been found at "$rcfile
#Add ersilia alias
if grep -q "alias ersilia" $rcfile; then
        echo "alias to ersilia already exists in " $rcfile
else
echo "alias ersilia=\"docker exec -it ersilia /bin/bash\"" >> $rcfile
echo "Alias added to ersilia"
alias ersilia="docker exec -it ersilia /bin/bash"
fi
if docker ps|grep -q ersilia; then
        echo "Ersilia is up and running!"
        docker ps|grep "ersilia"
        echo "To start using ersilia just type 'ersilia' in your shell."
else
        echo "Something went wrong!"
fi
