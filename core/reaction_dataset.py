import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import CountVectorizer

def load_and_preprocess_dataset(path="data/jee_reactions.csv"):
    df = pd.read_csv(path)

    # Combine reactants and products into one string
    df["reaction"] = df["reactants"] + " >> " + df["products"]
    
    X = df["reaction"]
    y = df["label"]

    vectorizer = CountVectorizer()
    X_vect = vectorizer.fit_transform(X)

    return train_test_split(X_vect, y, test_size=0.2, random_state=42), vectorizer
