input="bin/number1.txt"
while       read -r line; do
        number1="$line"
done < "$input"
echo $number1
