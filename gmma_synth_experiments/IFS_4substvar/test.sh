echo "Do you wish to overwrite output from previous run?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) overwrite; continue;;
        No ) exit 0;;
    esac 
done
echo "continued"
