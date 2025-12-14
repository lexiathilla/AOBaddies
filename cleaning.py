import pandas as pd

def remove_rows_with_asterisks(input_file, output_file, sheet_name="Sheet1", new_sheet_name="Cleaned"):
    # Load the Excel file
    df = pd.read_excel(input_file, sheet_name=sheet_name)
    print("file loaded")
    # Filter out rows where any cell contains '****'
    df_cleaned = df[~df.apply(lambda row: row.astype(str).str.contains(r"\*{4}").any(), axis=1)]
    
    # Save to a new sheet in the same workbook
    with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
        df_cleaned.to_excel(writer, sheet_name=new_sheet_name, index=False)

# Example usage
remove_rows_with_asterisks(r"C:\Users\Alexia\OneDrive - Imperial College London\AAYEAR4\APO1\GCs.xlsx",r"C:\Users\Alexia\OneDrive - Imperial College London\AAYEAR4\APO1\GCs_cleaned.xlsx")
