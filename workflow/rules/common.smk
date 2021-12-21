def annotate(tbl, annot, tbl_key, annot_key, annot_columns):
    annot_tbl = annot[annot_columns]
    annot_tbl = annot_tbl.drop_duplicates(subset=[annot_key])

    # merge tbl with annotation
    merged = tbl.merge(
        annot_tbl, left_on=tbl_key, right_on=annot_key, how="left", validate="1:1"
    )

    # reorder columns
    columns = (
        [tbl_key]
        + list(filter(annot_key.__ne__, annot_columns))
        + list(filter(tbl_key.__ne__, tbl.columns))
    )

    return merged[columns]
