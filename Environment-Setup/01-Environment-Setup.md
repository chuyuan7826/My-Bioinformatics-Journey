# ⚒️ This document records how to manage environment with `conda`
🚨 If you don't believe in me, copy and paste this page into the chat box of an AI you trust to see whether they are correct. I ⚠️**highly recommend**⚠️ you to do this, because I maintain this repo myself, so it is easy for me to make some mistakes.

## 1. What is a environment?

### 1.1 Problems encountered
In bioinformatics, command line tools, programming languages, and packages of **different versions** are very often to be used. If all the softwares are crammed into one folder, it will be quite difficult to figure out which version is desired to be invoked.

### 1.2 My definition
We need to create some separate folders to store the softwares of the specific versions. **One such folder is an environment**.

## 2. What is `conda`?

### 2.1 Core idea
`conda` is a **package manager** that can help create and swtich between these folders easily. In most cases, I use conda in the **terminal** as a **command line tool**.

### 2.2 Installation
Go to [Anaconda's official website](https://www.anaconda.com/download) and follow the instructions. I bet everyone can make it as long as you pay some patience. I recommend install **miniconda**, because it is lightweight and flexible.

## 3. Manage environment
Open the terminal, type `conda -V` and press enter. You will see return message that tells you the version of `conda` you installed. 

Type `conda env list` to see the environments in your disposal. You can see an asterisk (*) on the activated environment. `base` is created by default with the installation of miniconda.

If you don't see any asterisk, type `conda activate base` to activate `base` environment. Type `conda deactivate` to deactivate it.

Type `conda  create -n my_env_name` to create your own environment. `conda` will ask you to confirm the installation, type `y` (yes) to permit it. You can choose whatever name you like. Type `conda activate my_env_name` to activate the environment you just created.

Type `conda install python` (‼️ always run this command after you activated an environment) to install some packages (or you may want to call them softwares). `conda` will ask you to confirm the installation, type `y` (yes) to permit it.

Type `conda list` to see what packages are installed in your current environment. Type `conda list python` to filter for `python`.

Deactivate your environment (remember `conda deactivate`). Type `conda env remove -n my_env_name` to delete this environment. Don't worry, you can easily create environments by reuse the commands you've just learned.

## 4. Utilize AI
Ask an AI you trust the following questions:

* What are conda channels?
* How to export my conda  environment configuration?

Of course, you can ask anything you want.

## 5. 🎉 Explore working with conda yourself 🎉
Congratulations! You've stepped out further!

## 6. Quick Command Cheat Sheet
| Command | Purpose |
| :--- | :--- |
| `conda env list` | List all environments |
| `conda activate <name>` | Switch to an environment |
| `conda deactivate` | Leave current environment |
| `conda create -n <name>` | Create new environment |
| `conda install <package>` | Install a package |
| `conda list` | List installed packages |
