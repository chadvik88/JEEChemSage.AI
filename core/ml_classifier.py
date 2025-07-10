from sklearn.ensemble import RandomForestClassifier
from core.reaction_dataset import load_and_preprocess_dataset
import joblib
import os

MODEL_PATH = "core/model/rf_classifier.pkl"
VECT_PATH = "core/model/vectorizer.pkl"

def train_model():
    (X_train, X_test, y_train, y_test), vectorizer = load_and_preprocess_dataset()
    clf = RandomForestClassifier(n_estimators=150, random_state=42)
    clf.fit(X_train, y_train)

    acc = clf.score(X_test, y_test)
    print(f"Test Accuracy: {acc:.2f}")

    os.makedirs("core/model", exist_ok=True)
    joblib.dump(clf, MODEL_PATH)
    joblib.dump(vectorizer, VECT_PATH)
