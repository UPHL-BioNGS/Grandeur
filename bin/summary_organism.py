import pandas as pd

def predict_organism(summary_df):

    ##########################################
    # predicting organism                    #
    ##########################################

    print("Predicting organism")

    summary_df['predicted_organism'] = pd.NA

    if 'skani_organism' in summary_df.columns:
        summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna(summary_df['skani_organism'])

    if 'kraken2_top_organism' in summary_df.columns:
        summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna(summary_df['kraken2_top_organism'])

    if 'mash_screen_organism' in summary_df.columns:
        summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna(summary_df['mash_screen_organism'])

    if 'mash_dist_organism' in summary_df.columns:
        summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna(summary_df['mash_dist_organism'])

    summary_df['predicted_organism'] = summary_df['predicted_organism'].fillna("Unknown")

    summary_df['predicted_organism'] = summary_df['predicted_organism'].str.strip()

    # Adjusting shigellas
    if 'shigapass_predicted_serotype' in summary_df.columns:
        print("Refining E. coli and Shigella predictions using ShigaPass")
        
        def refine_shigella_ecoli(row):
            org = str(row.get('predicted_organism', ''))
            shiga_pred = str(row.get('shigapass_predicted_serotype', ''))
            
            if 'escherichia' in org.lower() or 'shigella' in org.lower():
                
                if pd.notna(shiga_pred) and shiga_pred.strip() != '' and shiga_pred.lower() != 'nan':
                    shiga_pred_clean = shiga_pred.strip().upper()
                    new_org = None
                    
                    if shiga_pred_clean.startswith('SS'):
                        new_org = 'Shigella sonnei'
                    elif shiga_pred_clean.startswith('SF'):
                        new_org = 'Shigella flexneri'
                    elif shiga_pred_clean.startswith('SB'):
                        new_org = 'Shigella boydii'
                    elif shiga_pred_clean.startswith('SD'):
                        new_org = 'Shigella dysenteriae'
                    elif 'NOT SHIGELLA' in shiga_pred_clean or 'EIEC' in shiga_pred_clean:
                        if 'shigella' in org.lower():
                            new_org = 'Escherichia coli'
                        else:
                            return org
                    
                    if new_org:
                        if '_' in org:
                            return new_org.replace(' ', '_')
                        else:
                            return new_org
                            
            return org

        summary_df['predicted_organism'] = summary_df.apply(refine_shigella_ecoli, axis=1)

    return summary_df