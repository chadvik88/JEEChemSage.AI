from core.ml_classifier import train_model

if __name__ == "__main__":
    train_model()

import pandas as pd

df = pd.read_csv("data/jee_reactions.csv")  # or whatever your dataset path is
print(df.head(10))
print("Unique labels:", df["label"].nunique())
print("Total samples:", len(df))
