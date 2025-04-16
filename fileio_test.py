
from fileio import load_globals

def main():
    filepath = "input_template.xlsx"  # Replace with full path if needed
    globals_data = load_globals(filepath)

    for key, value in globals_data.items():
        print(f"\n=== {key} ===")
        if isinstance(value, list):
            for item in value:
                print(item)
        else:
            print(value)

if __name__ == "__main__":
    main()
