#Importing Dependencies :- ///////////////////////////////////////////////////////////////////

import pandas as pd 
import numpy as np 
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler
import csv
import sys  
#Data Collection :- //////////////////////////////////////////////////////////////////////////


data = pd.read_csv("./CONTAMINATION_TRAINING_DATA.csv")


#Data preprocessing :- //////////////////////////////////////////////////////////////////////

x = data.drop(columns=['ID','LABEL'] , axis = 1)
y = data['LABEL']
x_train , x_test , y_train , y_test = train_test_split(x,y,stratify=y , test_size=0.2 , random_state=3)
scaler = StandardScaler()
scaler.fit(x_train.values)
x_train_std = scaler.transform(x_train.values)
x_test_std = scaler.transform(x_test.values)

#Model on the Board :- ///////////////////////////////////////////////////////////////////////
model = LogisticRegression()
model.fit(x_train_std , y_train)

#Import argumets :- //////////////////////////////////////////////////////////////////////////

ID = sys.argv[1]
Total_read_count = sys.argv[2]
Unmapped_count = sys.argv[3]
Mapped_Read_count = sys.argv[4]
Mapped_Read_coverage = sys.argv[5]
U_T = sys.argv[6]
Total_coverage = sys.argv[7]
Read_count = sys.argv[8]
microbe_count = sys.argv[9]
Norm = sys.argv[10]
Modest_microbe = sys.argv[11]
output_report = sys.argv[12]



#Parseing data :- ////////////////////////////////////////////////////////////////////////////

input_data = (Total_read_count,Unmapped_count,Mapped_Read_count,Mapped_Read_coverage,U_T,Total_coverage,Read_count,microbe_count,Norm,)
input_data_as_numpy_array = np.asarray(input_data)
input_data_reshaped = input_data_as_numpy_array.reshape(1,-1)
std_input = scaler.transform(input_data_reshaped)
prediction = model.predict(std_input)


#Writing report :- ///////////////////////////////////////////////////////////////////////////////
if prediction[0] == 0:
   status = "Yes"
else:
   status = "NO"

data = [
    ("ID", ID),
    ("Total Read Count", Total_read_count),
    ("Unmapped Read Count",Unmapped_count),
    ("Mapped Read Count" ,Mapped_Read_count ),
    ("Mapped Read Coverage" , Mapped_Read_coverage),
    ("U_T" , U_T),
    ("Total Coverage" ,Total_coverage ),
    ("Read Count of Microbe" ,Read_count ),
    ("Microbe Count" ,microbe_count ) , 
    ("Normalised read Count of Microbe" , Norm ),
    ("Contamination" , status ),
    ("Microbe with highest read count" , Modest_microbe)
    
]
columns = ["Feature", "Value"]

df = pd.DataFrame(data, columns=columns)

df.to_csv(output_report, sep='\t', quoting=csv.QUOTE_NONE, index=False, header=True)
