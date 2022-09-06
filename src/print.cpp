#include <string>
#include <vector>

#include "print.h"
#include "dist.h"

std::string GREEN(char c) { return "\033[32m" + std::string(1,c) + "\033[0m"; }
std::string GREEN(std::string str) { return "\033[32m" + str + "\033[0m"; }
std::string RED(char c) { return "\033[31m" + std::string(1,c) + "\033[0m"; }
std::string RED(std::string str) { return "\033[31m" + str + "\033[0m"; }
std::string BLUE(int i) { return "\033[34m" + std::to_string(i) + "\033[0m"; }
std::string BLUE(char c) { return "\033[34m" + std::string(1,c) + "\033[0m"; }
std::string BLUE(std::string str) { return "\033[34m" + str + "\033[0m"; }
std::string YELLOW(char c) { return "\033[33m" + std::string(1,c) + "\033[0m"; }
std::string YELLOW(std::string str) { return "\033[33m" + str + "\033[0m"; }
std::string PURPLE(char c) { return "\033[35m" + std::string(1,c) + "\033[0m"; }
std::string PURPLE(std::string str) { return "\033[35m" + str + "\033[0m"; }

void print_ptrs(std::vector< std::vector<int> > ptrs, 
        std::string alt_str, std::string ref_str) 
{

    // create array
    int alt_len = alt_str.size();
    int ref_len = ref_str.size();
    std::vector< std::vector<char> > ptr_str;
    for (int i = 0; i < alt_len*2; i++)
        ptr_str.push_back(std::vector<char>(ref_len*2, ' '));

    // set arrows
    for (int alt_idx = 0; alt_idx < alt_len; alt_idx++) {
        for(int ref_idx = 0; ref_idx < ref_len; ref_idx++) {
            // states
            if (ptrs[alt_idx][ref_idx] & PTR_DONE)
                ptr_str[alt_idx*2+1][ref_idx*2+1] = 'o';
            else if (ptrs[alt_idx][ref_idx] & PTR_NEXT)
                ptr_str[alt_idx*2+1][ref_idx*2+1] = '+';
            else
                ptr_str[alt_idx*2+1][ref_idx*2+1] = '*';

            // movements
            if (ptrs[alt_idx][ref_idx] & PTR_DIAG)
                ptr_str[alt_idx*2][ref_idx*2] = '\\';
            if (ptrs[alt_idx][ref_idx] & PTR_SUB)
                ptr_str[alt_idx*2][ref_idx*2] = 'X';
            if (ptrs[alt_idx][ref_idx] & PTR_LEFT)
                ptr_str[alt_idx*2+1][ref_idx*2] = '-';
            if (ptrs[alt_idx][ref_idx] & PTR_UP)
                ptr_str[alt_idx*2][ref_idx*2+1] = '|';
        }
    }

    // print array
    for (int i = -1; i < alt_len*2; i++) {
        for (int j = -1; j < ref_len*2; j++) {
            if (i < 0 && j < 0) {
                printf("\n  ");
            }
            else if (i < 0) {
                if (!(j%2)) printf("%c ", ref_str[j>>1]);
            } else if (j < 0) {
                if (!(i%2)) printf("%c ", alt_str[i>>1]);
                else printf("  ");
            } else {
                switch(ptr_str[i][j]) {
                    case 'o':
                    case '+':
                    case '*':
                        if (ptrs[i>>1][j>>1] & LEFT_PATH) {
                            if (ptrs[i>>1][j>>1] & RIGHT_PATH) { // both
                                printf("%s", GREEN(ptr_str[i][j]).data());
                            } else { // left
                                printf("%s", YELLOW(ptr_str[i][j]).data());
                            }
                        } else if (ptrs[i>>1][j>>1] & RIGHT_PATH) { // right
                            printf("%s", BLUE(ptr_str[i][j]).data());
                        } else { // none
                            printf("%c", ptr_str[i][j]);
                        }
                        break;
                    case 'X':
                    case '\\':
                        if (ptrs[i>>1][j>>1] & LEFT_PATH && 
                                (i>>1) > 0 && (j>>1) > 0 &&
                                ptrs[(i>>1)-1][(j>>1)-1] & LEFT_PATH) {
                            if (ptrs[i>>1][j>>1] & RIGHT_PATH &&
                                    ptrs[(i>>1)-1][(j>>1)-1] & RIGHT_PATH) { // both
                                printf("%s", GREEN(ptr_str[i][j]).data());
                            } else { // left
                                printf("%s", YELLOW(ptr_str[i][j]).data());
                            }
                        } else if (ptrs[i>>1][j>>1] & RIGHT_PATH &&
                                (i>>1) > 0 && (j>>1) > 0 &&
                                ptrs[(i>>1)-1][(j>>1)-1] & RIGHT_PATH) { // right
                            printf("%s", BLUE(ptr_str[i][j]).data());
                        } else { // none
                            printf("%c", ptr_str[i][j]);
                        }
                        break;
                    case '-':
                        if (ptrs[i>>1][j>>1] & LEFT_PATH && 
                                (j>>1) > 0 &&
                                ptrs[i>>1][(j>>1)-1] & LEFT_PATH) {
                            if (ptrs[i>>1][j>>1] & RIGHT_PATH &&
                                    ptrs[i>>1][(j>>1)-1] & RIGHT_PATH) { // both
                                printf("%s", GREEN('-').data());
                            } else { // left
                                printf("%s", YELLOW('-').data());
                            }
                        } else if (ptrs[i>>1][j>>1] & RIGHT_PATH && 
                                (j>>1) > 0 &&
                                ptrs[i>>1][(j>>1)-1] & RIGHT_PATH) { // right
                            printf("%s", BLUE('-').data());
                        } else { // none
                            printf("-");
                        }
                        break;
                    case '|':
                        if (ptrs[i>>1][j>>1] & LEFT_PATH && 
                                (i>>1) > 0 &&
                                ptrs[(i>>1)-1][j>>1] & LEFT_PATH) {
                            if (ptrs[i>>1][j>>1] & RIGHT_PATH &&
                                    ptrs[(i>>1)-1][j>>1] & RIGHT_PATH) { // both
                                printf("%s", GREEN('|').data());
                            } else { // left
                                printf("%s", YELLOW('|').data());
                            }
                        } else if (ptrs[i>>1][j>>1] & RIGHT_PATH && 
                                (i>>1) > 0 &&
                                ptrs[(i>>1)-1][j>>1] & RIGHT_PATH) { // right
                            printf("%s", BLUE('|').data());
                        } else { // none
                            printf("|");
                        }
                        break;
                    case ' ':
                            printf("%c", ptr_str[i][j]);
                        break;
                }
            }
        }
        printf("\n");
    }
}
