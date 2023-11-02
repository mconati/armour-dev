s = 1;
arr = [1:10:128, 128];
all_vertices = A.position_constraints(1, 121)
temp = size(all_vertices);
hold on
for i = 1:1:temp(2)
    link_vertices = all_vertices{i}
    k = boundary(link_vertices(:,1), link_vertices(:,2), link_vertices(:,3), s);
    h = trisurf(k,link_vertices(:,1),link_vertices(:,2),link_vertices(:,3),'FaceColor', 'r','FaceAlpha', 0.05,'EdgeColor', 'g','EdgeAlpha', 0.05);
end