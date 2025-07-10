import os

def get_mechanism_image_path(mechanism):
    if not mechanism:
        print("[DEBUG] No mechanism provided.")
        return None

    keyword = mechanism.split()[0].upper()  # E2 Elimination → E2
    filename = f"{keyword.lower()}.png"

    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'assets', 'mechanisms'))
    full_path = os.path.join(base_dir, filename)

    print(f"[DEBUG] Looking for mechanism image at: {full_path}")
    print(f"[DEBUG] Exists? {os.path.exists(full_path)}")

    return full_path if os.path.exists(full_path) else None
