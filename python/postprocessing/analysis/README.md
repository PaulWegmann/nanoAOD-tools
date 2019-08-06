Setting up a new release area:

    source /cvmfs/cms.cern.ch/cmsset_default.sh (you may want to add this to your .bashrc)
    cd /nfs/dust/cms/user/<youraccount>/
    cmsrel CMSSW_10_4_0
    cd CMSSW_10_4_0/src
    cmsenv
    git cms-init
    git clone https://github.com/adewit/nanoAOD-tools.git PhysicsTools/NanoAODTools
    cd PhysicsTools/NanoAODTools
    git checkout -b summer-project-2019 origin/summer-project-2019
    cd ../../
    scram b 


When you log in again, the release area will still be there but you have to repeat these steps:
    source /cvmfs/cms.cern.ch/cmsset_default.sh (unless you put this in your .bashrc to be executed by default when you log in)
    cd /nfs/dust/cms/user/<youraccount>/CMSSW_10_4_0/src/
    cmsenv


When you make changes to the code it's a good idea to push them to the repository. If you make an account on github, let me know your github username and I can add you to the repository as a collaborator.

For learning how to use git, this is a good starting point: https://www.codecademy.com/learn/learn-git

