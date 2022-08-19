#!/bin/bash

usersname=#FILLIN"First Last"
username=#FILLIN"githubusername"
email=#FILLIN"JHEDID@jhu.edu"
token=#FILLIN"token"
year=2022


#git config
git config --global user.name $usersname
git config --global user.email $email
git config --global commit.gpgsign false
git config --global core.editor "nano -w"

# Create qbbYEAR-answers directory
reponame=qbb${year}-answers
mkdir ~/${reponame}
cd ~/${reponame}
echo "# QBB${year} repository" > README.md
for word in *.bam *.bcf *.bed *.bg *.bw *.faa *.fastq *.fasta *.fa *.fai *.fna *.fq *.fra *.gff* *.gtf *.gz *.h5 *.sam *.seq *.vcf
do
  echo $word >> .gitignore
done

#git add directory as a repo
git init
git add .gitignore README.md
git commit -m "initial commit"
git branch -M main
git remote add origin https://github.com/${username}/${reponame}.git
curl -u ${username}:${token} https://api.github.com/user/repos -d "{\"name\": \"${reponame}\"}"
git push -u origin main --force

#Add collaborators
##method curl: https://docs.github.com/en/rest/collaborators/collaborators#add-a-repository-collaborator
for collaborator in andrew-bortvin cutsort dtaylo95 kweav msauria rmccoy7541 stephaniemyan
do
  curl \
    -X PUT \
    -H "Accept: application/vnd.github+json" \
    -H "Authorization: token ${token}" \
    "https://api.github.com/repos/${username}/${reponame}/collaborators/${collaborator}" \
    -d ''
done

# clone cmdb-quantbio repo
cd ~
git clone https://github.com/bxlab/cmdb-quantbio.git
