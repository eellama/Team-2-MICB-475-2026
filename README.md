
### Agenda for Feb 3 Meeting 
**Decide on a research question**
- HIV Dataset
	- Healthy vs HIV patient. In HIV infected persons, IL6 is chronically elevated. IL6 is known to play a key role in regulation of the microbiome. Thus, investigating differences in microbiome by IL6 levels can inform issues on comorbidities. 5 aims: 
	- Everything on project 2 must be done on our shared container. 
	- We MUST process the dataset prior to proposal submission (Feb 22)
		1. QIIME -> process HIV dataset
		- Filter out visits. 
			- They asked patients to come in 3 times, but only 2 and 3 did they collect fecal samples. These are technical replicates, we don’t want them in there: some patients are represented twice, others once.
			- FILTER FOR ONLY VISIT 2
		- IL6 vs HIV status (pos/neg) comparision
		2. alpha diversity metrics (run all of them)
			- Keep IL6 a continuous variable. 
			- IL6 vs some kind of alpha diversity. 
			- Hopefully, IL6 levels are positively correlated with diversity. 		    
			- pearsons
			- Not longitudinal, just correlation
			- 1 graph of Alpha diversity vs IL-6 in HIV+, another graph of alpha diversity vs IL-6 in HIV- (can be same graph dif colours)
		3. Beta diversity
			- Based on the outcome of Alpha diversity. Depending on which one has correlation this depends on which one we do. If it flattens out, maybe there’s a threshold where it doesn’t matter
			- Bin IL6 levels into High Medium Low
		4. Once binned… (WEEK 6 modules)
			1. If faith significant, do a core microbiome indicates taxa analysis. 
				1. Shows how many microbiomes unique to each other
				2. Shows taxa which are unique to a given condition
			2. If shannon significant, do a DESEQ analysis
				1. Means theres changes in abundance of microbes,
		5. If None of 4 significant
			1. Functional analysis using piecrust2
				1. Optional module. 
				2. This predicts the functional profile of our groups: what microbial pathways are upregulated and downregulated? 
			    
**what variables that may confound the data that we want to control for**
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
- team proposal
	- note on intro background
	- hypothesis/prediction
	- experimental aims/rationale
	- proposed approach
	- dataset overview 





