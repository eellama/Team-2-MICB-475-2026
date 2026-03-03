# Agenda for March 3 Meeting

### Pre-meeting notes
	- Follow up: How is the Alpha Diversity Analysis going? (Owen)
	- Follow up: How is the Beta Diversity Analysis going? (Florence)
	- Proposal feedback:
		- Noted on group chat that we intend to edit and resubmit to recover marks.
		- Any specific questions for Bessie?
		- Decide what edits are actually necessary vs unnecessary. 
			- Are we going to meet up sometime to discuss? Or individually focus on the sections we originally worked on?
	
### Meeting notes
	- 

### Action items
	- See Feb 24 Meeting for timeline. Kevin: Core microbiome analysis by X. Sarah: Indicator Species by X. 


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





