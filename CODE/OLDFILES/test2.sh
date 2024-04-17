line=0
input="CURRENT/number1.txt"
    while read -r line; do
         number="$line"
    done < "$input"
echo $number
