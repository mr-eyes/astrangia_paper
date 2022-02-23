tax_id=$1
echo $tax_id | taxonkit lineage | awk '$2!=""' | taxonkit reformat --format "{k},{p},{c},{o},{f},{g},{s},{t}" | cut -f 3