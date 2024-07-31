#OUT-OF-SAMPLE TIME ANALYSIS

##Load in, reference lists, and preprocessing

###Initialising reference objects for selecting relevant columns/indexes in functions
drops = ["AMP","SAM","TZP","CZO","CRO","CAZ","FEP","MEM","CIP","GEN","SXT","NIT","VAN" ]
to_drop = ["subject_id","micro_specimen_id"]

###Preprocessing
ref_df = pd.read_csv("urines_5t_1.csv")
ref_df['standard_age'] = ref_df['standard_age'].map(str)
dummy_df = pd.get_dummies(ref_df.drop(drops, axis=1))

###Initialising empty dictionaries and lists to populate
vari_overall = {}
c_values = {}
class_reps = {}

##Cross-validation: 2008-2010 training dataset

###Tested on 2008-2010 dataset
own_time_per("urines_5t_1.csv")
aucrocs_5t_1_1 = aucrocs
aucrocs_5t_1_1['df1'] = 1
aucrocs_5t_1_1['df2'] = 1
aucrocs_5t_1_1 = pd.DataFrame(aucrocs_5t_1_1)

###Tested on 2011-2013 dataset
across_time_per("urines_5t_1.csv","urines_5t_2.csv")
aucrocs_5t_1_2 = aucrocs
aucrocs_5t_1_2['df1'] = 1
aucrocs_5t_1_2['df2'] = 2
aucrocs_5t_1_2 = pd.DataFrame(aucrocs_5t_1_2)

###Tested on 2014-2016 dataset
across_time_per("urines_5t_1.csv","urines_5t_3.csv")
aucrocs_5t_1_3 = aucrocs
aucrocs_5t_1_3['df1'] = 1
aucrocs_5t_1_3['df2'] = 3
aucrocs_5t_1_3 = pd.DataFrame(aucrocs_5t_1_3)

###Tested on 2017-2019 dataset
across_time_per("urines_5t_1.csv","urines_5t_4.csv")
aucrocs_5t_1_4 = aucrocs
aucrocs_5t_1_4['df1'] = 1
aucrocs_5t_1_4['df2'] = 4
aucrocs_5t_1_4 = pd.DataFrame(aucrocs_5t_1_4)

##Cross-validation: 2011-2013 training dataset

###Tested on 2008-2010 dataset
across_time_per("urines_5t_2.csv","urines_5t_1.csv")
aucrocs_5t_2_1 = aucrocs
aucrocs_5t_2_1['df1'] = 2
aucrocs_5t_2_1['df2'] = 1
aucrocs_5t_2_1 = pd.DataFrame(aucrocs_5t_2_1)

###Tested on 2011-2013 dataset
own_time_per("urines_5t_2.csv")
aucrocs_5t_2_2 = aucrocs
aucrocs_5t_2_2['df1'] = 2
aucrocs_5t_2_2['df2'] = 2
aucrocs_5t_2_2 = pd.DataFrame(aucrocs_5t_2_2)

###Tested on 2014-2016 dataset
across_time_per("urines_5t_2.csv","urines_5t_3.csv")
aucrocs_5t_2_3 = aucrocs
aucrocs_5t_2_3['df1'] = 2
aucrocs_5t_2_3['df2'] = 3
aucrocs_5t_2_3 = pd.DataFrame(aucrocs_5t_2_3)

###Tested on 2017-2019 dataset
across_time_per("urines_5t_2.csv","urines_5t_4.csv")
aucrocs_5t_2_4 = aucrocs
aucrocs_5t_2_4['df1'] = 2
aucrocs_5t_2_4['df2'] = 4
aucrocs_5t_2_4 = pd.DataFrame(aucrocs_5t_2_4)

##Cross-validation: 2014-2016 training dataset

###Tested on 2008-2010 dataset
across_time_per("urines_5t_3.csv","urines_5t_1.csv")
aucrocs_5t_3_1 = aucrocs
aucrocs_5t_3_1['df1'] = 3
aucrocs_5t_3_1['df2'] = 1
aucrocs_5t_3_1 = pd.DataFrame(aucrocs_5t_3_1)

###Tested on 2011-2013 dataset
across_time_per("urines_5t_3.csv","urines_5t_2.csv")
aucrocs_5t_3_2 = aucrocs
aucrocs_5t_3_2['df1'] = 3
aucrocs_5t_3_2['df2'] = 2
aucrocs_5t_3_2 = pd.DataFrame(aucrocs_5t_3_2)

###Tested on 2014-2016 dataset
own_time_per("urines_5t_3.csv")
aucrocs_5t_3_3 = aucrocs
aucrocs_5t_3_3['df1'] = 3
aucrocs_5t_3_3['df2'] = 3
aucrocs_5t_3_3 = pd.DataFrame(aucrocs_5t_3_3)

###Tested on 2017-2019 dataset
across_time_per("urines_5t_3.csv","urines_5t_4.csv")
aucrocs_5t_3_4 = aucrocs
aucrocs_5t_3_4['df1'] = 3
aucrocs_5t_3_4['df2'] = 4
aucrocs_5t_3_4 = pd.DataFrame(aucrocs_5t_3_4)

##Cross-validation: 2017-2019 training dataset

###Tested on 2008-2010 dataset
across_time_per("urines_5t_4.csv","urines_5t_1.csv")
aucrocs_5t_4_1 = aucrocs
aucrocs_5t_4_1['df1'] = 4
aucrocs_5t_4_1['df2'] = 1
aucrocs_5t_4_1 = pd.DataFrame(aucrocs_5t_4_1)

###Tested on 2011-2013 dataset
across_time_per("urines_5t_4.csv","urines_5t_2.csv")
aucrocs_5t_4_2 = aucrocs
aucrocs_5t_4_2['df1'] = 4
aucrocs_5t_4_2['df2'] = 2
aucrocs_5t_4_2 = pd.DataFrame(aucrocs_5t_4_2)

###Tested on 2014-2016 dataset
across_time_per("urines_5t_4.csv","urines_5t_3.csv")
aucrocs_5t_4_3 = aucrocs
aucrocs_5t_4_3['df1'] = 4
aucrocs_5t_4_3['df2'] = 3
aucrocs_5t_4_3 = pd.DataFrame(aucrocs_5t_4_3)

###Tested on 2017-2019 dataset
own_time_per("urines_5t_4.csv")
aucrocs_5t_4_4 = aucrocs
aucrocs_5t_4_4['df1'] = 4
aucrocs_5t_4_4['df2'] = 4
aucrocs_5t_4_4 = pd.DataFrame(aucrocs_5t_4_4)

aucrocs_5t = pd.concat([aucrocs_5t_1_1, aucrocs_5t_1_2,
                     aucrocs_5t_1_3,aucrocs_5t_1_4,
                     aucrocs_5t_2_1,aucrocs_5t_2_2,
                     aucrocs_5t_2_3,aucrocs_5t_2_4,
                     aucrocs_5t_3_1,aucrocs_5t_3_2,
                     aucrocs_5t_3_3,aucrocs_5t_3_4,
                     aucrocs_5t_4_1,aucrocs_5t_4_2,
                     aucrocs_5t_4_3,aucrocs_5t_4_4],ignore_index=True)

aucrocs_5t.to_csv('aucrocs_5t.csv',index=False)

