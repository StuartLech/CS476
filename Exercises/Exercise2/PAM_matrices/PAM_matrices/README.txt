
README for PAM Matrix Generation

Username: slech
Operating System: macOS
Tools Used: Terminal (for SSH and SFTP)

---

1. File Transfer Process:
I transferred the Perl program `PAMCompMod.pl` and the input file `PAM1-mutprob.txt` from my local machine (Mac) to the home.cs.siue.edu server using SFTP. Below are the steps I followed:

- I used the `sftp` command in the Terminal to connect to the server:
  sftp slech@home.cs.siue.edu
- After successfully connecting to the server, I used the `put` command to upload the files from my Downloads folder:
  put /Users/startltech/Downloads/PAMCompMod.pl
  put /Users/startltech/Downloads/PAM1-mutprob.txt

2. Running the Perl Program:
Once the files were on the server, I connected to the server via SSH using the following command:
  ssh slech@home.cs.siue.edu
After logging in, I navigated to the directory where the files were stored and ran the Perl program using the command:
  perl PAMCompMod.pl

The program prompted me for a divergence value, and I input the following values to generate the PAM matrices:
  5, 10, 25, 50, 100, 250, 500, 1000, 2000, 3000

For each divergence value, the output was saved to a file named PAM<n>scores.txt (e.g., PAM5scores.txt for a divergence of 5).

3. Downloading the PAM Matrices:
After generating the matrices, I used SFTP to download the resulting files back to my local machine:
  sftp slech@home.cs.siue.edu
  get PAM5scores.txt
  get PAM10scores.txt
  ...
  (repeated for all 10 matrices)

4. Zipping the Files:
Once all the files were downloaded, I used the following command to zip them into a folder:
  zip -r PAM_matrices.zip PAM5scores.txt PAM10scores.txt ... README.txt

---
