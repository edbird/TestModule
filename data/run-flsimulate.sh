#!/bin/bash
#flsimulate -c simulation.conf --output-file flsimulate.brio
#/home/ecb/Falaise/build/BuildProducts/bin/flsimulate -c simulation.conf --output-file flsimulate.brio

function write_config
{

    #echo "NUMBER_OF_EVENTS=$NUMBER_OF_EVENTS"

    echo '#@description Simulation Test'
    echo '#@key_label "name"'
    echo '#@meta_label "type"'
    echo ''
    echo '[name="flsimulate" type="flsimulate::section"]'
    echo 'numberOfEvents : integer = '"$NUMBER_OF_EVENTS"
    echo ''
    echo '[name="flsimulate.variantService" type="flsimulate::section"]'
    echo 'settings : string[1] = "primary_events:generator=Se82.2nubb"'

}

function main
{

    ARG_C="$#"
    ARG_IX=1
    echo "ARG_C=$ARG_C"
    echo "ARGS: $@"

    # defaults
    NUMBER_OF_EVENTS=10000
    NUMBER_OF_EVENTS_FLAG="TRUE"

    RUN_FLAG="TRUE"

    if (( $ARG_C > 1 ))
    then
        while true

            ARG="${!ARG_IX}"
            echo "ARG_IX=$ARG_IX -> ARG=$ARG"

            if [ "$ARG" == "-n" ] || [ "$ARG" == "--number-of-events" ]
            then
                (( ARG_IX += 1 ))
                NUMBER_OF_EVENTS="${!ARG_IX}"
                (( ARG_IX += 1 ))
                echo "[ INFO ] : NUMBER_OF_EVENTS=$NUMBER_OF_EVENTS"

            elif [ "$ARG" == "--help" ] || [ "$ARG" == "-h" ]
            then
                (( ARG_IX += 1 ))
                echo "Usage $0 --number-of-events NUMBER_OF_EVENTS"
                echo ""
                echo "NUMBER_OF_EVENTS: Number of events to generate"

            else
                echo "[ WARNING ] : Unrecognized argument $ARG, ARG_IX=$ARG_IX"
                echo "[ INFO ] : Run $0 --help for help"
                (( ARG_IX += 1 ))

            fi

        do

            # Standard continue / break sequence
            if (( ARG_IX <= ARG_C ))
            then
                continue
            else
                break
            fi

        done
    else
        echo "Usage: $0 -n NUMBER_OF_EVENTS"
        RUN_FLAG="FALSE"
    fi

    
    if [ "$NUMBER_OF_EVENTS_FLAG" == "TRUE" ] && [ "$RUN_FLAG" == "TRUE" ]
    then
        # Edit flsimualate.conf
        #cat simulation.conf.backup
        write_config > ./flsimulate.conf
        #cat simulation.conf
        #diff simulation.conf.backup simulation.conf

        # Run job
        /home/ecb/Falaise/build/BuildProducts/bin/flsimulate -c flsimulate.conf --output-file flsimulate.brio
    
    else
        echo "[ ERROR ] : Invalid NUMBER_OF_EVENTS"

    fi
}

main "$@"
