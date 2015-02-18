# 
# chr12:3109870-3110051
# chr14:19415634-19419744
# chr2:98661401-98668291
# chr9:2998303-3039411
# chr9:35305007-35305778

#Snps in the region
chromosome.region.1 <- subset(wholeLhSnpList, wholeLhSnpList$loc > 3109870 & wholeLhSnpList$loc < 3110051)

chromosome.region.2 <- subset(wholeLhSnpList, wholeLhSnpList$loc > 19415634 & wholeLhSnpList$loc < 19419744)

chromosome.region.3 <- subset(wholeLhSnpList, wholeLhSnpList$loc > 98661401 & wholeLhSnpList$loc < 98668291)

chromosome.region.4 <- subset(wholeLhSnpList, wholeLhSnpList$loc > 2998303 & wholeLhSnpList$loc < 3039411)

chromosome.region.5 <- subset(wholeLhSnpList, wholeLhSnpList$loc > 35305007 & wholeLhSnpList$loc < 35305778)

