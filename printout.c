#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

int main() {
	int fd;
	char c;
	
	int i=0;
	fd = open("SPEC110215", O_RDONLY);
	while ( (read(fd, &c, 1) > 0)) {
		//printf("byte %d, %d\n",i, c);
		i++;
	}
	printf("%d bytes total\n",i);
	close(fd);
}
