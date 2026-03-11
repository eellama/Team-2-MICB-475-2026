# Agenda for March 10 Meeting
### Pre-meeting notes
- Proposal revision(?)
- Follow up: Beta Diversity Analysis results (Florence)
- Follow up: Core Microbiome Analysis results (Kevin)
- Follow up: Indicator Species Analysis results (Sarah)
- Follow up: DESeq Analysis results (Elena)
- Clarify: Do we still want to include healthy controls for downstream analyses?
- Clarify: Ritu suggested we do PICRUSt2 analysis on the entire microbiome instead of just ASVs identified in aim 3 of the proposal. Do we want to do this?
  
### Meeting notes
- HIV_all affected IL-6, but effects are minimal. Thus, we will focus on HIV+ only
- Decide if we want to glom
- Core mcirobiome and ISA results align
- Consider increasing threshold for Core Microbiome

### Action items
- crosscheck indicator species to lit review
- Check if 2-3 samples with higher IL-6 levels identified in HIV_all & HIVpos beta analyses have significantly higher IL-6 than the rest of the sample (Florence)
- Titrate which increased threshold works for core microbiome
- Glom to genus level for core microbiome, keep at 0.7, see what we will find. Consider (but don't do) anything to prevalence think about it
- Glom before creating volcano plot for DESeq

# Agenda for March 3 Meeting

### Pre-meeting notes
- Follow up: How is the Alpha Diversity Analysis going? (Owen)
- Follow up: How is the Beta Diversity Analysis going? (Florence)
- Proposal feedback:

  -> Noted on group chat that we intend to edit and resubmit to recover marks.

- Identify any specific questions for Bessie.

  - In last week's meeting, we discussed binning IL-6 levels of healthy controls as well, but this was not flagged in the proposal. Should we do our analyses based on 3 experimental groups (not binning IL-6 of healthy controls) or 4 (binning IL-6 of healthy controls too)?

    -> Beta analysis done with 4 experimental groups
- Decide what edits are actually necessary vs unnecessary.

  -> Are we going to meet up sometime to discuss? Or should we individually focus on the sections we originally worked on?
	
### Meeting notes
- Alpha diversity update:
  - Results from Faith's PD and Shannon diversity
  - If Faith's PD significant -> Core microbiome analysis & indicator species analysis
  - If Shannon significant -> DESeq differencial abundance analysis

- Beta diversity update:
  - PCoA/clustering patterns

- May be better to focus on Beta diversity. Owen's quick analysis shows no significant difference in alpha diversity.

### Action items
- See Feb 24 Meeting for overall timeline.
  - Data processing: Internal deadline by Sunday
  - Kevin: Core microbiome analysis by March 10. 
  - Sarah: Indicator Species by March 10. 
  - Elena: Deseq by March 10. 
- Group: Review TA feedback and prepare proposal revisions
  - Assign sections & set deadline for revision?
- Proposal revision due Sunday internal. Use the google docs SUGGESTING FEATURE to show that we are addressin ghte comments


# Agenda for Feb 24 Meeting

### Notes
- Rarified phyloseq should ONLY be for diversity metrics.
- all others use unrarified.
- New idea: bin IL-6 High and Low amongst Healthy and HIV+ persons (4 unique categories). Then, in beta diversity, see which is more correlated with the other? If its IL-6, then look at that. If it's HIV, then look at that.
- If the TA feedback doesnt dock marks for the proposal don't change anything. Don't make unnecessary edits, edit according to TA feedback. 

### Clarify data processing roles, next step
- Kevin
- Owen
- Florence
- Elena
-  Sarah (Absent due to illness)

### Internal deadlines
- **26th Feb**, Thursday: Phyloseq object created by Elena
- **2nd Feb**, Monday: Beta diversity analysis by Florence, Alpha diversity analysis by Owen

### Action Items
- Do AGENDA AT LEAST 12-24 HOURS BEFORE THE MEETING!!!
- Workflow from today
	- Phyloseq object creation by Elena due by Thursday 26th -> distribute to group google drive
	- Beta Diversity by Florence by March 2
	- Alpha Diversity by Owen by March 2
- AFTER MEETING MARCH 3
	- Core Microbiome by Kevin after Diversity metrics are out
	- Indicator Species by Sarah after Diversity metrics are out
	- Differential gene abundance by X
	- Functional analysis by Owen
- Overall workflow
	1. Week of March 2: Aim 1, Revise proposal, address all TA comments to get points back
	2. Week of March 9: Aim 2
	3. Week of March 16: Aim 3
	4. Week of March 23-30: presentation and drafting 


# Agenda for Feb 10 Meeting
**Follow up: Data processing subteam**
- How is that going
- ........ the train left the station 

**Follow up: Lit review/writing team**
- how is that going

### Action Items
- clarify what statistical tests to conduct for each type of analysis
- clarify how to best represent results of each analysis
- what needs to be cited (how cite programs used)
- decide exactly what to filter for
- clarify how to consider effects of treatment
- what beta diversity functions to use
- do we need a project timeline in the proposal?
- *internal* (is there a better way to format this readme i.e. separate pre- and post-meeting notes)
- Lit review people: look at if there has been past literature discussing HIV, IL-6, and impacts on microbiome
- Lit review: Short paragraph to describe the original dataset

# Agenda/Notes for Feb 3 Meeting 
**Decide on a research question:**
- HIV Dataset
	- Healthy vs HIV patient.
 	- In HIV infected persons, IL6 is chronically elevated.
  	- IL6 is known to play a key role in regulation of the microbiome.
  	- Thus, investigating differences in microbiome by IL6 levels can inform issues on comorbidities. 5 aims: 
 - Notes:
 	- Everything on project 2 must be done on our shared container. 
	- We MUST process the dataset prior to proposal submission (Feb 22)
   	- What variables that may confound the data that we want to control for?

**5 Aims**
1. QIIME -> process HIV dataset
	- Filter out visits. FOR ONLY VISIT 2
		- They asked patients to come in 3 times, but only collected fecal samples on day 2+3. 
		- These are technical replicates; some patients are represented twice, others once.
	- IL6 vs HIV status (pos/neg) comparision
2. alpha diversity metrics (run all of them)
	- Keep IL6 a continuous variable. 
	- IL6 vs some kind of alpha diversity. 
	- Hopefully, IL6 levels are positively correlated with diversity. 		    
	- pearsons
	- Not longitudinal, just correlation
	- 1 graph of Alpha diversity vs IL-6 in HIV+, another graph of alpha diversity vs IL-6 in HIV- (can be same graph dif colours)
3. Beta diversity
	- Based on the outcome of Alpha diversity. Pick which one to do depending on which one has correlation. 
	- If it flattens out, maybe there’s a threshold where it doesn’t matter
	- Bin IL6 levels: High Medium Low
4. Once binned… (WEEK 6 modules)
	- If faith significant, do a core microbiome indicates taxa analysis. 
		- Shows how many microbiomes unique to each other
		- Shows taxa which are unique to a given condition
	- If shannon significant, do a DESEQ analysis
		- Means theres changes in abundance of microbes,
5. If None of 4. is significant (can still do even if 4 is sig)
	- Functional analysis using piecrust2
		- Optional module. 
		- This predicts the functional profile of our groups: what microbial pathways are upregulated and downregulated?

Experimental conditions:
- neg control
- pos w/o treatment
- pos w/ treatment

**Meeting notes**
- Evelyn will send group server credentials tomorrow morning. 

**Pre-meeting notes**
- Do we want to use Zotero on our google docs? Yes please
- I made a document for the group lit review. When we are doing our lit review, please take notes on the team proposal document ‘lit review tab’ and upload a pdf of the document on the lit review folder. 
- Please change the name of the paper on the pdf name
- Please group it in the appropriate sub-category folder
- Team proposal is due FEB 22 Midnight (18 days from today)

**Action Items…**
- Begin literature review, using the documentation system described above
- Download dataset from server and process data by X 
- familiarize with the planned analyses
- begin proposal
	- Proposal expectations for next week
- contact Bessie before each step (email/canvas)
- WHAT ARE WE BRINGING NEXT WEEK?-
- Lit Review, Kevin and Sarah
	- HIV reviews (Sarah) 
	- IL6 in context of the immune system (Kevin)
	- IL6 and microbiome impact  (Florence)
	- Research gap 
- Data Processing:
	- Demultiplexing
	- Denoising/Cluttering
	- Excel IL6 HIV significance 
- Team Proposal:
	- notes on intro background
	- hypothesis/prediction
	- experimental aims/rationale
	- proposed approach
	- dataset overview 





