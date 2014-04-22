def make_horizontal_links(width, height, mod):
	links = []
	for y in range(0, height):
		for x in range(0, width, 2):
			link_id = y * width + x + mod
			if (x + mod + 1 < width):
				links.append((link_id, link_id + 1))
			
	return links
	
def make_vertical_links(width, height, mod):
	links = []
	for x in range(0, width):
		for y in range(0, height, 2):
			link_id = (y + mod) * width + x
			if (y + mod + 1 < height):
				links.append((link_id, link_id + width))
			
	return links
	
def make_forward_diagonal_links(width, height, mod):
	links = []
	for x in range(0, width):
		for y in range(0, height, 2):
			link_id = (y + mod) * width + x
			if (x + 1 < width and y + mod + 1 < height):
				links.append((link_id, link_id + 1 + width))
			
	return links
	
def make_backward_diagonal_links(width, height, mod):
	links = []
	for x in range(0, width):
		for y in range(0, height, 2):
			link_id = (y + mod) * width + x + 1
			if (x + 1 < width and y + mod + 1 < height):
				links.append((link_id, link_id + width - 1))
			
	return links
	
def make_horizontal_bend_links(width, height, x):
	links = []
	for y in range(0, height):
		link_id = y * width + x
		links.append((link_id, link_id + 2))
			
	return links
	
def make_vertical_bend_links(width, y):
	links = []
	for x in range(0, width):
		link_id = y * width + x
		links.append((link_id, link_id + 2 * width))
			
	return links
	
def make_links(width, height):
	links = []
	links.append(make_horizontal_links(width, height, 0))
	links.append(make_horizontal_links(width, height, 1))
	# print links
	
	links.append(make_vertical_links(width, height, 0))
	links.append(make_vertical_links(width, height, 1))
	
	links.append(make_forward_diagonal_links(width, height, 0))
	links.append(make_forward_diagonal_links(width, height, 1))
	
	links.append(make_backward_diagonal_links(width, height, 0))
	links.append(make_backward_diagonal_links(width, height, 1))
	
	for x in range(0, width - 2):
		links.append(make_horizontal_bend_links(width, height, x))
	
	for y in range(0, height - 2):
		links.append(make_vertical_bend_links(width, y))
		
	return links
		
if __name__ == "__main__":
	width = 6
	height = 4
	links = make_links(width, height)
	print len(links)
	
	