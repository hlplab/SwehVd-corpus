#########################################################################################
# ARPABET VOWEL ANALYZER                                                                #
#                                                                                       #
# DESCRIPTION:                                                                          #
# This script (modeled on Mietta Lennes' collect_formant_data_from_files.praat          #
# available at http://www.helsinki.fi/~lennes/praat-scripts/ and distributed under the  # 
# GNU General Public License, copyright 4/7/2003) is designed to be run on a set of     #
# soundfiles and P2FA-generated TextGrids. It extracts duration (in ms), timestamps     #
# (in s), and F1/F2/F3 from all intervals containing Arpabet vowels, and extracts       #
# labels for corresponding word and preceding & following phones, along with any notes  #
# present in the TextGrid. Last, the script outputs date, settings,   			#
# OS, and Praat version to the results file.The script can also be constrained to a     #
# user-defined set of words using the "targets" option. Before writing the results      #
# file, words are lowercased and Arpabet is converted to Unicode IPA.                   #
#                                                                                       #
# NOTES:                                                                                #
# 1) Remember to add a / or \ (depending on the OS) to the end of all paths. 2)         #
# Soundfiles and TextGrids must have identical names in order for the script to match   #
# them up. 3) To use the "targets" option, make sure your TextGrid includes a word tier #
# and the "use word tier" option is selected, then create a tab-delimited text file     #
# containing a list of the words (case-sensitive) you'd like to extract, separated by   #
# newlines, with "word" as the column header. 4) Because this script converts Arpabet   #
# to Unicode IPA, make sure Praat's text writing preferences are set to UTF-8 before    #
# running the script. If you intend to view the results file in Excel, in Windows you   #
# will have to open it from within Excel and import it specifying UTF-8 encoding in     #
# order for the characters to display correctly. On OSX, you will first have to convert #
# the file to UTF-16 Little Endian using a text editor, then import it into Excel.      #
#                                                                                       #
# CHANGELOG: 
# 10/10/22: Added possibility to measure at 5 time points                               #
# 11/12/20: Added token label to each analyzed vowel.					#
# 04/18/15: Fixed bug in code preceding/following phone code, fixed error in            #
#           instructions.                                                               #
# 02/08/14: Rewrote code for results file generation, added pitch extraction option,    #
#           added metadata output to the script (analyst, version, settings, etc.). 	#
# 11/06/13: Corrected ARPABET-IPA conversion (thanks Daniel).                           #
# 06/26/13: Fixed bug introduced by a recent version of Praat.                          #
# 04/27/13: Fixed error in Arpabet -> IPA conversion.                                   #
# 03/10/13: Reordered formant options, fixed preceding/following phone extraction bug.  #
# 02/22/13: Added option to append data to an existing results file.                    #
# 02/21/13: Added counter for number of vowels analyzed if using targets option, and    #
#           the option to only analyze stressed vowels.                                 #
# 02/14/13: Fixed bugs involving running the script over multiple files, and extracting #
#           surrounding phonological environment. Also simplified the extraction of     #
#           sounds from longsounds.                                                     #
# 01/19/13: Added ability to select the formant measurement points.                     #
# 01/01/13: Release version.                                                            #
#                                                                                       #
# This modified script distributed under the GNU General Public License v3 or higher,   #
# copyright 1/2013, John Riebold (riebold@uw.edu).                                      #
# 											#
# NOTES BY AP:										#	
# The converting to IPA-command is commented out.					#
# 											#
#########################################################################################

# PROMPT THE USER FOR THE LOCATION OF THE INPUT/OUTPUT FILES, FORMANT SETTINGS, ETC.
form Arpabet Vowel Analyzer
	sentence Results_file vowel_formants_swe.txt
	optionmenu Use_targets_file 2
		option yes
		option no
	sentence Targets_file targets.txt
	comment Tiers:
	sentence Phone_tier phone
	sentence Word_tier word
	optionmenu Use_notes_tier: 2
		option yes
		option no
	sentence Notes_tier notes
	comment Options:
	optionmenu Analyze_unstressed_vowels: 2
		option yes
		option no
	comment Formant settings:
	optionmenu Measurement_points: 6
		option Midpoint
		option 35%/50%/65%
		option 30%/50%/70%
		option 25%/50%/75%
		option 20%/50%/80%
		option 20%/35%/50%/65%/80%
	positive Maximum_formant_(Hz) 5500
	integer Number_of_formants 5
	comment Pitch Settings
	optionmenu Extract_pitch: 2
		option yes
		option no
	integer left_Pitch_range_(Hz) 75
	integer right_Pitch_range_(Hz) 500
endform

# SOUNDFILE DIRECTORY
soundfile_directory$ = chooseDirectory$: "Choose the directory for the sound files"
textgrid_directory$ = chooseDirectory$: "Choose the directory for the textgrids"
soundfile_directory$ = soundfile_directory$ + "/"
textgrid_directory$ = textgrid_directory$ + "/"

# SET ADDITIONAL FORMANT OPTIONS, CHANGE IF NECESSARY
preemphasis_from = 50
window_length = 0.025
time_step = 0.01

# DEFINE EMPTY VARIABLE IN CASE LABEL EMPTY/NOT PRESENT IN TEXTGRID
notes_label$ = ""

# DEFINE DUMMY COUNTER VARIABLES FOR END-OF-SCRIPT REPORT
sound_count = 0
vowel_count = 0
target_vowel_count = 0

# GET TIME AND OS
rundate$ = date$ ()
if windows = 1
    os$ = "Windows"
elsif macintosh = 1
    os$ = "OSX"
elsif unix = 1
	os$ = "Linux"
endif
version$ = "'praatVersion'"
version$ = replace_regex$ ("'version$'", "(\d)(\d)(\d{2,2})", "\1.\2.\3", 0)

# PRESUMABLY ADDED THIS LINE SO THAT PRAAT SHOULD STORE OUTPUT IN CORRECT DIR, BUT IT'S NOT WORKING
#results_file$ = "../../data/" + results_file$

# INITIALIZE RESULTS FILE
if fileReadable (results_file$)
	beginPause ("Warning")
		comment ("The file 'results_file$' already exists.")
	results_choice = endPause ("Append", "Overwrite", 1)
	if results_choice = 2
		filedelete 'results_file$'
		call InitializeResultsFile
	endif
else
	call InitializeResultsFile
endif

# OPEN TARGETS FILE
if use_targets_file = 1
	Read Table from tab-separated file... 'targets_file$'
	targets$ = selected$ ("Table", 1)
endif

# CREATE LIST OF SOUNDFILES IN DIRECTORY
Create Strings as file list... list 'soundfile_directory$'*.wav
numberoffiles = Get number of strings

# GO THROUGH EACH SOUND FILE
for ifile to numberoffiles
	select Strings list
	filename$ = Get string... ifile

	# OPEN SOUNDFILE FROM LIST
	Open long sound file... 'soundfile_directory$''filename$'
	soundfile$ = selected$ ("LongSound", 1)

	# INCREMENT SOUND COUNT
	sound_count = sound_count + 1

	# OPEN TEXTGRID OF SAME NAME
	gridfile$ = "'textgrid_directory$''soundfile$'.TextGrid"
	if fileReadable (gridfile$)
		Read from file... 'gridfile$'

		# FIND TIER NUMBER FOR PHONE AND WORD TIERS
		call GetTier 'phone_tier$' phone_tier
		call GetTier 'word_tier$' word_tier
		intervals = Get number of intervals... phone_tier

		# EXTRACT ANNOTATED PORTION OF SOUNDFILE
		gridstart = Get start time
		gridend = Get end time
		select LongSound 'soundfile$'
		Extract part... gridstart gridend yes

		# REMOVE LONGSOUND
		select LongSound 'soundfile$'
		Remove

		# EXTRACT FORMANT AND PITCH OBJECTS
		select Sound 'soundfile$'
		To Formant (burg)... time_step number_of_formants maximum_formant window_length preemphasis_from
		if extract_pitch = 1
			select Sound 'soundfile$'
			To Pitch... 0 left_Pitch_range right_Pitch_range
		endif

		# PASS THROUGH EACH INTERVAL IN SELECTED TIER AND GET LABEL
		for interval to intervals
			select TextGrid 'soundfile$'
			phone_label$ = Get label of interval... phone_tier interval

			# CHECK IF INTERVAL CONTAINS ARPABET VOWEL, IF SO, ANALYZE IT
			if analyze_unstressed_vowels = 1
				if phone_label$ = "AA1" or phone_label$ = "AH1" or phone_label$ = "EE1" or phone_label$ = "EH1" or phone_label$ = "II1" or phone_label$ = "IH1" or phone_label$ = "OO1" or phone_label$ = "OH1" or phone_label$ = "UU1" or phone_label$ = "UH1" or phone_label$ = "YY1" or phone_label$ = "YH1" or phone_label$ = "AE1" or phone_label$ = "{:" or phone_label$ = "{" or phone_label$ = "OA1" or phone_label$ = "OAH1" or phone_label$ = "OE1" or phone_label$ = "OEH1" or phone_label$ = "9:" or phone_label$ = "9"
				call AnalyzeVowel
				endif
			elsif analyze_unstressed_vowels = 2
				if phone_label$ = "AA1" or phone_label$ = "AH1" or phone_label$ = "EE1" or phone_label$ = "EH1" or phone_label$ = "II1" or phone_label$ = "IH1" or phone_label$ = "OO1" or phone_label$ = "OH1" or phone_label$ = "UU1" or phone_label$ = "UH1" or phone_label$ = "YY1" or phone_label$ = "YH1" or phone_label$ = "AE1" or phone_label$ = "{:" or phone_label$ = "{" or phone_label$ = "OA1" or phone_label$ = "OAH1" or phone_label$ = "OE1" or phone_label$ = "OEH1" or phone_label$ = "9:" or phone_label$ = "9"
				call AnalyzeVowel
				endif
			endif
		endfor

		# REMOVE TEXTGRID OBJECT FROM OBJECT LIST
		select TextGrid 'soundfile$'
		Remove
	endif

	# REMOVE TEMPORARY OBJECTS AND CONTINUE WITH NEXT FILE
	select Sound 'soundfile$'
	plus Formant 'soundfile$'
	if extract_pitch = 1
		plus Pitch 'soundfile$'
	endif
	Remove
endfor

# REMOVE REST OF OBJECTS AND FINISH
select Strings list
if use_targets_file = 1
	plus Table 'targets$'
endif
Remove

# PRINT A REPORT
echo Done. Analyzed 'target_vowel_count' of 'vowel_count' vowels in 'sound_count' file(s).

# PROCEDURE TO ANALYZE VOWELS
procedure AnalyzeVowel

	# INCREMENT VOWEL COUNT
	vowel_count = vowel_count + 1

	# GET START AND END TIMES, CALCULATE DURATION, ETC.
	start = Get starting point... phone_tier interval
	end = Get end point... phone_tier interval
	duration = (end-start)
	duration_ms = duration*1000
	midpoint = (start+end)/2
	
	# DETERMINE WHICH POINTS TO MEASURE
	if measurement_points = 2
		onset = start+(duration*0.35)
		offset = start+(duration*0.65)
	elsif measurement_points = 3
		onset = start+(duration*0.30)
		offset = start+(duration*0.70)
	elsif measurement_points = 4
		onset = start+(duration*0.25)
		offset = start+(duration*0.75)
	elsif measurement_points = 5
		onset = start+(duration*0.20)
		offset = start+(duration*0.80)
	elsif measurement_points = 6
		onset = start+(duration*0.20)
		point2 = start+(duration*0.35)
		point4 = start+(duration*0.65)
		offset = start+(duration*0.80)
	endif

	# GET FORMANT VALUES AT INTERVAL(S)
	select Formant 'soundfile$'
	f1_2 = Get value at time... 1 midpoint Hertz Linear
	f2_2 = Get value at time... 2 midpoint Hertz Linear
	f3_2 = Get value at time... 3 midpoint Hertz Linear
	if measurement_points = 6
		f1_1 = Get value at time... 1 onset Hertz Linear
		f2_1 = Get value at time... 2 onset Hertz Linear
		f3_1 = Get value at time... 3 onset Hertz Linear
		f1_3 = Get value at time... 1 point2 Hertz Linear
		f2_3 = Get value at time... 2 point2 Hertz Linear
		f3_3 = Get value at time... 3 point2 Hertz Linear
		f1_4 = Get value at time... 1 point4 Hertz Linear
		f2_4 = Get value at time... 2 point4 Hertz Linear
		f3_4 = Get value at time... 3 point4 Hertz Linear
		f1_5 = Get value at time... 1 offset Hertz Linear
		f2_5 = Get value at time... 2 offset Hertz Linear
		f3_5 = Get value at time... 3 offset Hertz Linear
	elsif measurement_points > 1 and measurement_points < 6
		f1_1 = Get value at time... 1 onset Hertz Linear
		f2_1 = Get value at time... 2 onset Hertz Linear
		f3_1 = Get value at time... 3 onset Hertz Linear
		f1_3 = Get value at time... 1 offset Hertz Linear
		f2_3 = Get value at time... 2 offset Hertz Linear
		f3_3 = Get value at time... 3 offset Hertz Linear
	endif
	
	# EXTRACT PITCH AT INTERVAL(S)
	if extract_pitch = 1
		select Pitch 'soundfile$'
		f0 = Get mean... midpoint-0.02 midpoint+0.02 Hertz
	endif

	# GET WORD VOWEL IS FROM
	select TextGrid 'soundfile$'
	word = Get interval at time... word_tier midpoint
	word_label$ = Get label of interval... word_tier word

	# GET PRECEDING AND FOLLOWING ENVIRONMENTS, SKIPPING SPACES
	preceding_label$ = Get label of interval... phone_tier (interval-1)
	if preceding_label$ = "sp" or preceding_label$ = "sil"
		if interval-2 >= 1
			preceding_label$ = Get label of interval... phone_tier (interval-2)
		elsif interval-2 < 1
			preceding_label$ = ""
		endif
	endif
	if interval <= intervals-1
		following_label$ = Get label of interval... phone_tier (interval+1)
		if following_label$ = "sp" or following_label$ = "sil"
			if interval+2 <= intervals
				following_label$ = Get label of interval... phone_tier (interval+2)
			elsif interval+2 > intervals
				following_label$ = ""
			endif
		endif
	elsif interval+1 > intervals
		following_label$ = ""
	endif

	# INSERT TOKEN LABEL IN ORDER TO IDENTIFY INDIVIDUAL TOKENS IN DATA
	if phone_label$ = "AA1" or phone_label$ = "AH1" or phone_label$ = "EE1" or phone_label$ = "EH1" or phone_label$ = "II1" or phone_label$ = "IH1" or phone_label$ = "OO1" or phone_label$ = "OH1" or phone_label$ = "UU1" or phone_label$ = "UH1" or phone_label$ = "YY1" or phone_label$ = "YH1" or phone_label$ = "AE1" or phone_label$ = "{:" or phone_label$ = "{" or phone_label$ = "OA1" or phone_label$ = "OAH1" or phone_label$ = "OE1" or phone_label$ = "OEH1" or phone_label$ = "9:" or phone_label$ = "9"
		for i from 1 to vowel_count
			token_label$ = phone_label$ + "_" + string$ (i)
		endfor
	endif
	
	# GET CONTENTS OF NOTES TIER
	if use_notes_tier = 1
		call GetTier 'notes_tier$' notes_tier
		note = Get interval at time... notes_tier midpoint
		notes_label$ = Get label of interval... notes_tier note
	endif

	# CONVERT WORDS TO LOWERCASE, ARPABET TO UNICODE IPA
	word_label$ = replace_regex$ (word_label$, "[A-Z]", "\L&", 0)
	call ConvertText phone_label$
	call ConvertText preceding_label$
	call ConvertText following_label$

	# CREATE RESULTS LINE
	resultsline_begin$ = "'soundfile$'	'word_label$'	'phone_label$'	'token_label$'	'preceding_label$'	'following_label$'	'start'	'end'	'duration'	"
	if use_notes_tier = 1
		resultsline_end$ = "'notes_label$'	'rundate$'	Max formant: 'maximum_formant' Hz, Number of formants: 'number_of_formants', Window length: 'window_length' s	'version$'	'os$''newline$'"
	else
		resultsline_end$ = "'rundate$'	Max formant: 'maximum_formant' Hz, Number of formants: 'number_of_formants', Window length: 'window_length' s	'version$'	'os$''newline$'"
	endif
	if measurement_points = 1
		resultsline_middle$ = "'f1_2'	'f2_2'	'f3_2'	"
	elsif measurement_points = 6
		resultsline_middle$ = "'f1_1'	'f1_3'	'f1_2'	'f1_4'	'f1_5'	'f2_1'	'f2_3'	'f2_2'	'f2_4'	'f2_5'	'f3_1'	'f3_3'	'f3_2'	'f3_4'	'f3_5'	"
	elsif measurement_points > 1 and measurement_points < 6
		resultsline_middle$ = "'f1_1'	'f1_2'	'f1_3'	'f2_1'	'f2_2'	'f2_3'	'f3_1'	'f3_2'	'f3_3'	"		
	endif
	if extract_pitch = 1
		resultsline_middle$ = "'f0'	" + resultsline_middle$
	endif
	resultsline$ = resultsline_begin$ + resultsline_middle$ + resultsline_end$

	# OUTPUT TO RESULTS FILE
	if use_targets_file = 1
		select Table 'targets$'
		match = Search column... word 'word_label$'
		if match
			target_vowel_count = target_vowel_count + 1
			fileappend "'results_file$'" 'resultsline$'
		endif
	else
		target_vowel_count = target_vowel_count + 1
		fileappend "'results_file$'" 'resultsline$'
	endif
endproc

# PROCEDURE TO INITIALIZE RESULTS FILE
procedure InitializeResultsFile
	header_begin$ = "Filename	Word	Vowel	Token	Preceding Phone	Following Phone	Begin Time (s)	End Time (s)	Duration (s)	"
	if use_notes_tier = 1
		header_end$ = "Notes	Date	Settings	Praat Version	OS'newline$'"
	else
		header_end$ = "Date	Settings	Praat Version	OS'newline$'"
	endif
	if measurement_points = 1
		header_middle$ = "F1_50	F2_50	F3_50	"
	elsif measurement_points = 2
		header_middle$ = "F1_35	F1_50	F1_65	F2_35	F2_50	F2_65	F3_35	F3_50	F3_65	"
	elsif measurement_points = 3
		header_middle$ = "F1_30	F1_50	F1_70	F2_30	F2_50	F2_70	F3_30	F3_50	F3_70	"
	elsif measurement_points = 4
		header_middle$ = "F1_25	F1_50	F1_75	F2_25	F2_50	F2_75	F3_25	F3_50	F3_75	"
	elsif measurement_points = 5
		header_middle$ = "F1_20	F1_50	F1_80	F2_20	F2_50	F2_80	F3_20	F3_50	F3_80	"
	elsif measurement_points = 6
		header_middle$ = "F1_20	F1_35	F1_50	F1_65	F1_80	F2_20	F2_35	F2_50	F2_65	F2_80	F3_20	F3_35	F3_50	F3_65	F3_80	"
	endif
	if extract_pitch = 1
		header_middle$ = "F0.mean	" + header_middle$
	endif
	header$ = header_begin$ + header_middle$ + header_end$
	fileappend "'results_file$'" 'header$'
endproc

# PROCEDURE TO FIND NUMBER OF TIER WITH GIVEN LABEL
procedure GetTier name$ variable$
	numberOfTiers = Get number of tiers
	itier = 1
	repeat
		tier$ = Get tier name... itier
		itier = itier + 1
	until tier$ = name$ or itier > numberOfTiers
	if tier$ <> name$
		'variable$' = 0
	else
		'variable$' = itier - 1
	endif
	if 'variable$' = 0
		exit The tier 'name$' is missing from the file 'soundfile$'!
	endif
endproc

# PROCEDURE TO CONVERT ARPABET TO UNICODE IPA
procedure ConvertText arplabel$
	'arplabel$' = replace_regex$ ('arplabel$', "[A-Z]", "\L&", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "ao\d", "ɔ", 0) 
	#'arplabel$' = replace_regex$ ('arplabel$', "aa\d", "ɑ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "iy\d", "i", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "uw\d", "u", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "eh\d", "ɛ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "ih\d", "ɪ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "uh\d", "ʊ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "ah[12]", "ʌ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "ah0", "ə", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "ax\d", "ə", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "ae\d", "æ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "ey\d", "e", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "ay\d", "aj", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "ow\d", "o", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "aw\d", "aw", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "oy\d", "ɔj", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "er\d", "ɝ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "axr\d", "ɚ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "ch", "ʧ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "jh", "ʤ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "th", "θ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "dh", "ð", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "sh", "ʃ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "zh", "ʒ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "hh", "h", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "em", "m̩", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "en$", "n̩", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "^ng", "ŋ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "eng", "ŋ̩", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "el", "ɫ̩", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "r", "ɹ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "dx", "ɾ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "nx", "ɾ̃", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "q", "ʔ", 0)
	#'arplabel$' = replace_regex$ ('arplabel$', "y", "j", 0)
endproc