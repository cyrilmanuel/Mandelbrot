#include <stdio.h>
#include <stdlib.h>

int main()
{
    int n,m;
    printf("Entrer le nombre de thread :");
    scanf("%d",&n);
    printf("Entrer le nombre de bloc :");
    scanf("%d",&m);
    if(n>0 && m>0 && n<=m)
    {
        char str1[10], str2[10];
        char progCmdline[100];
        sprintf(str1, " %d", n);
        sprintf(str2, " %d", m);
        strcpy(progCmdline,"testSequentiel.exe ");
        strcat (progCmdline,str1);
        strcat (progCmdline,str2);
        system(progCmdline);
    }
    else
    {
        exit(0);
    }

}
