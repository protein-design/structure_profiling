'''
asynchronous running
Decoupled from bio_omics (bioomics.common.async_concurr),
used by the af_msa pipeline.
'''
import asyncio

async def run_worker(func, data_queue, result_queue):
    """define task and wrap it into queue"""
    while True:
        data = await data_queue.get()

        if data is None:
            data_queue.task_done()
            return
        
        # Send result to result queue
        try:
            result = await func(data)
            await result_queue.put(result) 
        finally:
            data_queue.task_done()


async def process_results(result_queue, total_items, post_func=None):
    """
    Consumer the processes result in living time
    """
    processed = 0
    while processed < total_items:
        result = await result_queue.get()
        try:
            # Process immediately
            if post_func is not None:
                post_func(result)
            processed += 1
        finally:
            result_queue.task_done()
        
async def main(task_funcs:list, data:list, post_func=None):
    '''
    example: asyncio.run(async_main(funcs, data))
    '''
    data_queue = asyncio.Queue()
    result_queue = asyncio.Queue()
    
    # create input queue
    for item in data:
        await data_queue.put(item)

    # distribute worker with task
    workers = []
    for func in task_funcs:
        task = asyncio.create_task(
            run_worker(func, data_queue, result_queue)
        )
        workers.append(task)

    # Start result processor
    result_processor = asyncio.create_task(
        process_results(result_queue, len(data), post_func)
    )
    
    # Complete task: wait till all queued tasks done
    await data_queue.join()

    # complete and shutdown a worker if no data
    for _ in workers:
        await data_queue.put(None)

    # wait for all workers complete
    await asyncio.gather(*workers)

    # Wait for all results to be processed
    await result_processor