# Day 5 Practical

This practical will cover:

1. Generating whole reference genome alignments
2. Visualising dot plots
3. Investigate structural rearrangements

We will investigate sex chromosome inversions identified in the fourspine stickleback by **[Liu et al. 2025](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011465)**

## 00. Prepare work folder for day 5
```
mkdir day5
cd day1
conda activate sexchr
```

## 01. Local installation of visualisation software

* **[JBrowse](https://jbrowse.org/jb2/)** - The next-generation genome browser

Download JBrowse, a standalone application, to your local machine, chose the apropriate version for your system platform,
double click the installation and follow the on screen instructions, double click the JBrowse App Icon

```
https://jbrowse.org/jb2/download/
```


## 02. Identify sex chromosome assemblies from INSDC
Search **[NCBI Genomes](https://www.ncbi.nlm.nih.gov/genome/)** for the entry "Apeltes quadracus"  
Select assembly **[GCA_048569185.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_048569185.1/)**  
Download the **[X](https://www.ncbi.nlm.nih.gov/nuccore/CM109091.1?report=fasta)** and the **[Y](https://www.ncbi.nlm.nih.gov/nuccore/CM109092.1?report=fasta)** chromosome fasta files to your local machine  

         
## 03. Align X and Y chromosome
Align the X and Y chromosome sequence online with **[DGenies](https://dgenies.toulouse.inra.fr)**  
Click on the "RUN" tab  


<img width="1293" height="831" alt="Screenshot 2025-09-29 at 17 12 08" src="https://github.com/user-attachments/assets/bbc23361-b8be-4c76-b17b-817a79c2941d" />

  
Select the fasta files you downloaded previously and fill in the required information.  
Click **Submit**  

Your status will change from  

  
<img width="817" height="362" alt="Screenshot 2025-09-29 at 17 17 08" src="https://github.com/user-attachments/assets/3a5cab80-3197-4b2d-9ef4-7472a95e4dc1" />  


To  

<img width="785" height="460" alt="Screenshot 2025-09-29 at 17 18 53" src="https://github.com/user-attachments/assets/0432353d-228d-41eb-8bf1-07703a0b0e0f" />

in a couple of minutes  

## 04. Inspect the results online

<img width="787" height="746" alt="Screenshot 2025-09-30 at 09 51 35" src="https://github.com/user-attachments/assets/1b957ff2-a4a2-4428-b35c-3cd015cde012" />  

## 05. Export the results
Export the results in PAF format to your local machine.  


## 06. Load the results into Jbrowse
Open Jbrowse  
Click **Open Sequence file(s)**  

<img width="622" height="634" alt="Screenshot 2025-09-30 at 10 01 54" src="https://github.com/user-attachments/assets/30c49768-3805-4381-8288-fe74f3c24020" />  

Click **ADD ANOTHER ASSEMBLY**  
  
<img width="617" height="691" alt="Screenshot 2025-09-30 at 10 03 25" src="https://github.com/user-attachments/assets/dee3ba26-7fad-4c33-b76c-d55de07755a6" />

Click **SUBMIT**  

Select  **Dotplot View -> LAUNCH VIEW**  
The following error message will appear  
<img width="1348" height="422" alt="Screenshot 2025-09-30 at 10 22 43" src="https://github.com/user-attachments/assets/96b6b184-6b0e-4a87-a10b-b49ed979e25f" />

Choose the following settings  



Click **LAUNCH**  
and explore the plot

Click **ADD** , choose **Linear Synteny View**  

<img width="1332" height="240" alt="Screenshot 2025-09-30 at 10 25 39" src="https://github.com/user-attachments/assets/f21e387c-7661-4d15-8c16-bae7d913c74d" />

Click **LAUNCH**  

You can now explore this plot  

<img width="1375" height="334" alt="Screenshot 2025-09-30 at 10 29 15" src="https://github.com/user-attachments/assets/da17496b-2731-4e18-8c35-482236cc63c0" />


